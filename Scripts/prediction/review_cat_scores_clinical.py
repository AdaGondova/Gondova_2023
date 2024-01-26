## requires ohbm env
import pandas as pd 
import numpy as np 
from sklearn.svm import NuSVR, NuSVC
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import LeaveOneOut, KFold, StratifiedKFold
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.metrics import mean_absolute_error
import pingouin as pg
import matplotlib.pyplot as plt 
import seaborn as sns
import pickle
import random
import pwlf
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, balanced_accuracy_score, mean_absolute_error, recall_score
import os


#=== GLOBAL INFO SETUP ===#
with open(r"../../DerivedData/cohorts_subjects_list.pickle", "rb") as input_file:
        cohorts = pickle.load(input_file)
        
# in the end outcomes are the scores, however here I am setting up the baselines predicting prematurity status etc
outcomes = ['Cognitive','Language','Motor', 'group']
#outcome='group_cat'
n_folds = 46 # to recreate train/test proportions 
it_num = 100 # for null distributions
num_rep = 10 # for repeated k-fold



#========================================================= READ IN DATA ============================================================#

df = pd.read_csv('../for_risk_and_env_factors_new_minAD.csv', index_col=0)

#============================================== categorise Prematurity prediction ==================================================#
df['group'] = 'FT'
df.loc[df['GA_birth'] < 37, 'group'] = 'PT'

df['group_cat'] = 0
df.loc[df['GA_birth'] < 37, 'group_cat'] = 1

for score in ['Cognitive', 'Language', 'Motor']: 
	df[score+'_cat'] = 0
	df.loc[df[score] <= 85, score+'_cat'] = 1
outcomes_cols = [score+'_cat' for score in outcomes]

'''inputs = {
    
    'GA birth' : [['GA_birth'], [0]], # 0 do not correct age, 1 do correct
    'PMA scan' : [['PMA_scan'], [0]],
    'FA' : [['FA'], [0]],
    'FA corr' : [['FA'], [1]],
    'GA birth + PMA scan' : [['GA_birth', 'PMA_scan'], [0,0]],
    'GA birth + PMA scan + FA' : [['GA_birth', 'PMA_scan', 'FA'], [0,0,0]],
    'GA birth + PMA scan + FA corrected' : [['GA_birth', 'PMA_scan', 'FA'], [0,0,1]],
    'Segmentation (52)' : [ FA_cols[:-1], np.zeros_like( FA_cols[:-1], dtype=int)],
    'Segmentation (52) corrected' : [ FA_cols[:-1], np.ones_like( FA_cols[:-1], dtype=int)],
}'''

inputs = {
	'clinical' : [['corr_FA', 'sex', 'GA_birth', 'IMD'], [0,0,0,0]],
	'no_age' : [['corr_FA', 'sex', 'IMD'], [0,0,0]]
}



#=============================================FUNCTIONS ===========================================================#
def get_null_distribution(X, y, PMA_scan, corr_idx, n_folds, mode, it=1000):
    #this is basically repeated K-fold where y is shuffled!
 
    if mode == 'cat':
        res = {'auc' : [], 'acc' : [], 'spec' : [], 'sens' : []}
        for i in range(it):
            y_shuff = np.random.permutation(y)
            y_true, y_pred = run_Kfold(X=X, y=y_shuff, PMA_scan=PMA_scan, 
                               corr_idx=for_corr_index, n_folds=n_folds, mode='cat')
            auc, acc, spec, sens = evaluate(y_true=y_true, y_pred=y_pred, mode='cat')
            res['auc'].append(auc)
            res['acc'].append(acc)
            res['spec'].append(spec)
            res['sens'].append(sens)
            
    elif mode == 'cont':
        res = {'rho' : [], 'pval' : [], 'mae' : []}
        for i in range(it):
            y_shuff = np.random.permutation(y)
            y_true, y_pred = run_Kfold(X=X, y=y_shuff, PMA_scan=PMA_scan, 
                               corr_idx=for_corr_index, n_folds=n_folds, mode='cont')
            rho, pval, mae = evaluate(y_true=y_true, y_pred=y_pred, mode='cont')
            res['rho'].append(rho)
            res['pval'].append(pval)
            res['mae'].append(mae)
    return res
    

def get_repeated_Kfold(X, y, PMA_scan, corr_idx, n_folds, mode, num_reps=100):   
    if mode == 'cat':
        res = {'auc' : [], 'acc' : [], 'spec' : [], 'sens' : []}
        for i in range(num_reps):
            X,y,PMA_scan = unison_shuffled_copies(a=X, b=y, c=PMA_scan)
            y_true, y_pred = run_Kfold(X=X, y=y, PMA_scan=PMA_scan, 
                               corr_idx=for_corr_index, n_folds=n_folds, mode='cat')
            auc, acc, spec, sens = evaluate(y_true=y_true, y_pred=y_pred, mode='cat')
            res['auc'].append(auc)
            res['acc'].append(acc)
            res['spec'].append(spec)
            res['sens'].append(sens)
            
    elif mode == 'cont':
        res = {'rho' : [], 'pval' : [], 'mae' : []}
        for i in range(num_reps):
            X,y,PMA_scan = unison_shuffled_copies(a=X, b=y, c=PMA_scan)
            y_true, y_pred = run_Kfold(X=X, y=y, PMA_scan=PMA_scan, 
                               corr_idx=for_corr_index, n_folds=n_folds, mode='cont')
            rho, pval, mae = evaluate(y_true=y_true, y_pred=y_pred, mode='cont')
            res['rho'].append(rho)
            res['pval'].append(pval)
            res['mae'].append(mae)
    return res
    
    
def unison_shuffled_copies(a, b, c):
    assert len(a) == len(b) == len(c)
    p = np.random.permutation(len(a))
    return a[p], b[p], c[p]
        
        
def run_Kfold(X, y, PMA_scan, corr_idx, n_folds, mode):
    
    loo = KFold(n_splits=n_folds)
    
    y_true = []
    y_pred = []
    
    for train_index, test_index in loo.split(X):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]
        PMA_train, PMA_test = PMA_scan[train_index], PMA_scan[test_index]
        
        ## process inputs 
        X_train, X_test = _pre_process(X_train=X_train, X_test=X_test, 
                                       PMA_train=PMA_train, PMA_test=PMA_test, 
                                       corr_idx=corr_idx, inflection=36)
        
        y_out = run_model(X_train=X_train, X_test=X_test, 
                           y_train=y_train, y_test=y_test, 
                           mode='cat')
        for el in range(len(y_test)):
                y_true.append(y_test[el])
                y_pred.append(y_out[el]) 
    return y_true, y_pred
        
def evaluate(y_true, y_pred, mode):
    
    if mode == 'cat':
        ## auc, balanced acc, specificity, sensitivity 
        return roc_auc_score(y_true, y_pred), balanced_accuracy_score(y_true, y_pred), recall_score(y_true, y_pred, pos_label=0),recall_score(y_true, y_pred, pos_label=1) 
    if mode == 'cont':
        return pg.corr(y_true, y_pred)['r'][0], pg.corr(y_true, y_pred)['p-val'][0], mean_absolute_error(y_true, y_pred)
        

def run_model(X_train, X_test, y_train, y_test, mode): 
    ## mode is 'cat' or 'cont'
    if mode == 'cat':
        clf = SVC( C=9, kernel='linear')
        clf.fit(X_train, y_train)
        y_out = clf.predict(X_test)
        
        y_out[y_out>0.5]=1
        y_out[y_out<=0.5]=0
    
    if mode == 'cont':
        clf = NuSVR(C=9, kernel='linear')
        clf.fit(X_train, y_train)
        y_out = clf.predict(X_test)
    
    return y_out

def _pre_process(X_train, X_test, PMA_train, PMA_test, corr_idx, inflection=36):
    
    # impute median
    X_train, X_test = _impute_median(X_train=X_train, X_test=X_test)
    # correct age
    if len(corr_idx) > 0 :
        X_train, X_test = _correct_age(X_train=X_train, X_test=X_test, 
                                      PMA_train=PMA_train, PMA_test=PMA_test, corr_idx=corr_idx, inflection=36)
        
    X_train, X_test = _scaling(X_train= X_train, X_test = X_test)
    
    return X_train, X_test
    

def get_columns_to_correct(inDict, key):
    in_cols, for_corr_index = inDict[key][0], np.argwhere(np.array(inDict[key][1]) == 1).ravel()
    if len(for_corr_index) > 0:
        print('Correcting inputs: {}'.format(np.array(in_cols)[for_corr_index]))   
    return in_cols, for_corr_index

def _impute_median(X_train, X_test):
    
    for col in range(len(X_train[0])):
        
        md = np.nanmedian(X_train[:,col])
        
        X_train[:,col][np.isnan(X_train[:,col])] = md
        X_test[:,col][np.isnan(X_test[:,col])] = md
    
    return X_train, X_test

def _scaling(X_train, X_test):
    
    scaler = MinMaxScaler()
    scaler.fit(X_train)
    
    X_train = scaler.transform(X_train)
    X_test = scaler.transform(X_test)
    
    return X_train, X_test

def _correct_age(X_train, X_test, PMA_train, PMA_test, corr_idx, inflection=36):
    x0 = np.array([min(PMA_train), inflection, max(PMA_train)])
    
    for idx in corr_idx:
        myPWLF = pwlf.PiecewiseLinFit(PMA_train, X_train[:,idx])
        myPWLF.fit_with_breaks(x0)
        
        ## correct train 
        yHat_train = myPWLF.predict(PMA_train)
        res_train = X_train[:,idx] - yHat_train
        
        ## correct test 
        yHat_test = myPWLF.predict(PMA_test)
        res_test = X_test[:,idx] - yHat_test
        
        X_train[:,idx] = res_train
        X_test[:,idx] = res_test
    return X_train, X_test   


# ================================================================== RUN PREDICTION ==================================================================#
o_file = '../../Results/predictions/review_cat_outcomes_on_IMD_results_minAD.txt'
if os.path.exists(o_file):
    os.remove(o_file)



### stupidest ever but I cannot survive another loop 
# ===================== COGNITIVE =======#
outcome='Cognitive_cat'
for cohort in ['A']:

	with open(o_file, 'a') as the_file:
		the_file.write('Cohort {}\n'.format(cohort))
		the_file.write('{}\n'.format(outcome))
		the_file.write('Inputs, AUC mean (std), ACC mean (std), Spec mean (std), Sens mean (std), Permuted AUC (p-val)\n')
	
	sub_df = df.copy()
	sub_df = sub_df.sample(frac=1).reset_index(drop=True)
	
	print('Cohort {}'.format(cohort))
	print('{}'.format(outcome))
	print('Inputs, AUC mean (std), ACC mean (std), Spec mean (std), Sens mean (std), Permuted AUC (p-val)')

	for key in inputs.keys():
    
		in_cols, for_corr_index = get_columns_to_correct(inDict=inputs, key=key)
    
    		### get inputs
		X = sub_df[in_cols].values
		y = sub_df[outcome].values
		PMA_scan = sub_df['PMA_scan'].values
    
    		### RUN & EVALUATE 
		res = get_repeated_Kfold(X=X, y=y, PMA_scan=PMA_scan, 
                             corr_idx=for_corr_index, n_folds=n_folds, 
                             mode='cat', num_reps=num_rep)
		null = get_null_distribution(X=X, y=y, PMA_scan=PMA_scan, 
                             corr_idx=for_corr_index, n_folds=n_folds, 
                             mode='cat', it=it_num)

		auc_p = np.sum(null['acc']>=res['acc'][0]) / len(null['acc'])

		out_message = '{}, {:.3f}({:.3f}),{:.3f}({:.3f}),{:.3f}({:.3f}),{:.3f}({:.3f}),{:.3f}({:.6f})\n'.format(
                        key,
                        np.mean(res['auc']), np.std(res['auc']),
                        np.mean(res['acc']), np.std(res['acc']),
                        np.mean(res['spec']), np.std(res['spec']),
                        np.mean(res['sens']), np.std(res['sens']),
                        res['auc'][0],auc_p
    			)
		print(out_message)
		with open(o_file, 'a') as the_file:
    			the_file.write(out_message)



# ===================== Language =======#
outcome='Language_cat'
for cohort in ['A']:

	with open(o_file, 'a') as the_file:
		the_file.write('Cohort {}\n'.format(cohort))
		the_file.write('{}\n'.format(outcome))
		the_file.write('Inputs, AUC mean (std), ACC mean (std), Spec mean (std), Sens mean (std), Permuted AUC (p-val)\n')

	sub_df = df.copy()
	sub_df = sub_df.sample(frac=1).reset_index(drop=True)

	
	print('Cohort {}'.format(cohort))
	print('{}'.format(outcome))
	print('Inputs, AUC mean (std), ACC mean (std), Spec mean (std), Sens mean (std), Permuted AUC (p-val)')

	for key in inputs.keys():
    
		in_cols, for_corr_index = get_columns_to_correct(inDict=inputs, key=key)
    
    		### get inputs
		X = sub_df[in_cols].values
		y = sub_df[outcome].values
		PMA_scan = sub_df['PMA_scan'].values
    
    		### RUN & EVALUATE 
		res = get_repeated_Kfold(X=X, y=y, PMA_scan=PMA_scan, 
                             corr_idx=for_corr_index, n_folds=n_folds, 
                             mode='cat', num_reps=num_rep)
		null = get_null_distribution(X=X, y=y, PMA_scan=PMA_scan, 
                             corr_idx=for_corr_index, n_folds=n_folds, 
                             mode='cat', it=it_num)

		auc_p = np.sum(null['acc']>=res['acc'][0]) / len(null['acc'])

		out_message = '{}, {:.3f}({:.3f}),{:.3f}({:.3f}),{:.3f}({:.3f}),{:.3f}({:.3f}),{:.3f}({:.6f})\n'.format(
                        key,
                        np.mean(res['auc']), np.std(res['auc']),
                        np.mean(res['acc']), np.std(res['acc']),
                        np.mean(res['spec']), np.std(res['spec']),
                        np.mean(res['sens']), np.std(res['sens']),
                        res['auc'][0],auc_p
    			)
		print(out_message)
		with open(o_file, 'a') as the_file:
    			the_file.write(out_message)
		

# ===================== Motor =======#
outcome='Motor_cat'
for cohort in ['A']:

	with open(o_file, 'a') as the_file:
		the_file.write('Cohort {}\n'.format(cohort))
		the_file.write('{}\n'.format(outcome))
		the_file.write('Inputs, AUC mean (std), ACC mean (std), Spec mean (std), Sens mean (std), Permuted AUC (p-val)\n')


	sub_df = df.copy()
	sub_df = sub_df.sample(frac=1).reset_index(drop=True)
	
	print('Cohort {}'.format(cohort))
	print('{}'.format(outcome))
	print('Inputs, AUC mean (std), ACC mean (std), Spec mean (std), Sens mean (std), Permuted AUC (p-val)')

	for key in inputs.keys():
    
		in_cols, for_corr_index = get_columns_to_correct(inDict=inputs, key=key)
    
    		### get inputs
		X = sub_df[in_cols].values
		y = sub_df[outcome].values
		PMA_scan = sub_df['PMA_scan'].values
    
    		### RUN & EVALUATE 
		res = get_repeated_Kfold(X=X, y=y, PMA_scan=PMA_scan, 
                             corr_idx=for_corr_index, n_folds=n_folds, 
                             mode='cat', num_reps=num_rep)
		null = get_null_distribution(X=X, y=y, PMA_scan=PMA_scan, 
                             corr_idx=for_corr_index, n_folds=n_folds, 
                             mode='cat', it=it_num)

		auc_p = np.sum(null['acc']>=res['acc'][0]) / len(null['acc'])

		out_message = '{}, {:.3f}({:.3f}),{:.3f}({:.3f}),{:.3f}({:.3f}),{:.3f}({:.3f}),{:.3f}({:.6f})\n'.format(
                        key,
                        np.mean(res['auc']), np.std(res['auc']),
                        np.mean(res['acc']), np.std(res['acc']),
                        np.mean(res['spec']), np.std(res['spec']),
                        np.mean(res['sens']), np.std(res['sens']),
                        res['auc'][0],auc_p
    			)
		print(out_message)
		with open(o_file, 'a') as the_file:
    			the_file.write(out_message)
	
# ===================== COGNITIVE =======#
outcome='group_cat'
for cohort in ['A']:

	with open(o_file, 'a') as the_file:
		the_file.write('Cohort {}\n'.format(cohort))
		the_file.write('{}\n'.format(outcome))
		the_file.write('Inputs, AUC mean (std), ACC mean (std), Spec mean (std), Sens mean (std), Permuted AUC (p-val)\n')
	
	sub_df = df.copy()
	sub_df = sub_df.sample(frac=1).reset_index(drop=True)
	
	print('Cohort {}'.format(cohort))
	print('{}'.format(outcome))
	print('Inputs, AUC mean (std), ACC mean (std), Spec mean (std), Sens mean (std), Permuted AUC (p-val)')

	for key in inputs.keys():
    
		in_cols, for_corr_index = get_columns_to_correct(inDict=inputs, key=key)
    
    		### get inputs
		X = sub_df[in_cols].values
		y = sub_df[outcome].values
		PMA_scan = sub_df['PMA_scan'].values
    
    		### RUN & EVALUATE 
		res = get_repeated_Kfold(X=X, y=y, PMA_scan=PMA_scan, 
                             corr_idx=for_corr_index, n_folds=n_folds, 
                             mode='cat', num_reps=num_rep)
		null = get_null_distribution(X=X, y=y, PMA_scan=PMA_scan, 
                             corr_idx=for_corr_index, n_folds=n_folds, 
                             mode='cat', it=it_num)

		auc_p = np.sum(null['acc']>=res['acc'][0]) / len(null['acc'])

		out_message = '{}, {:.3f}({:.3f}),{:.3f}({:.3f}),{:.3f}({:.3f}),{:.3f}({:.3f}),{:.3f}({:.6f})\n'.format(
                        key,
                        np.mean(res['auc']), np.std(res['auc']),
                        np.mean(res['acc']), np.std(res['acc']),
                        np.mean(res['spec']), np.std(res['spec']),
                        np.mean(res['sens']), np.std(res['sens']),
                        res['auc'][0],auc_p
    			)
		print(out_message)
		with open(o_file, 'a') as the_file:
    			the_file.write(out_message)

