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
from sklearn.metrics import roc_auc_score, balanced_accuracy_score, mean_absolute_error, recall_score, r2_score
import os


#=== GLOBAL INFO SETUP ===#
with open(r"../../DerivedData/cohorts_subjects_list.pickle", "rb") as input_file:
        cohorts = pickle.load(input_file)
        
# in the end outcomes are the scores, however here I am setting up the baselines predicting prematurity status etc
outcomes = ['Cognitive Score','Language Score','Motor Score']
#outcome='group_cat'
n_folds = 46 # to recreate train/test proportions 
it_num = 100 # for null distributions
num_rep = 10 # for repeated k-fold



#========================================================= READ IN DATA ============================================================#
# ages
df = pd.read_csv('../../DerivedData/cohortA_subjects_clinical.csv', index_col=0)

# global FA
glob = pd.read_csv('../../DerivedData/extracted_metrics/global_cortical_FA_metrics_median_minAD.csv', index_col=0)

metrics = ['FA']
hemispheres = ['left', 'right']
for metric in metrics:
    for i, row in glob.iterrows():
        glob.loc[i,metric] = np.mean([row['left_{}'.format(metric)], 
                                   row['right_{}'.format(metric)]])
df = pd.merge(df, glob[['subject_id', 'FA']], on=['subject_id'])  

# median FA regions 
diff = pd.read_csv('../../DerivedData/extracted_metrics/neonat_segmentation_FA_median_minAD.csv', index_col=0)


FA_cols = [col for col in diff.columns if 'FA' in col]
FA_cols.extend(['subject_id'])

#new_df = df.copy()
df = pd.merge(df, diff[FA_cols], on=['subject_id'])

#============================================== Inputs ==================================================#

inputs = {
    
    'GA birth' : [['GA_birth'], [0]], # 0 do not correct age, 1 do correct
    'PMA scan' : [['PMA_scan'], [0]],
    #'FA' : [['FA'], [0]],
    'FA corr' : [['FA'], [1]],
    'GA birth + PMA scan' : [['GA_birth', 'PMA_scan'], [0,0]],
    #'GA birth + PMA scan + FA' : [['GA_birth', 'PMA_scan', 'FA'], [0,0,0]],
    'GA birth + PMA scan + FA corrected' : [['GA_birth', 'PMA_scan', 'FA'], [0,0,1]],
    #'Segmentation (52)' : [ FA_cols[:-1], np.zeros_like( FA_cols[:-1], dtype=int)],
    'Segmentation (52) corrected' : [ FA_cols[:-1], np.ones_like( FA_cols[:-1], dtype=int)],
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
        res = {'rho' : [], 'pval' : [], 'mae' : [], 'r2' : []}
        for i in range(it):
            y_shuff = np.random.permutation(y)
            y_true, y_pred = run_Kfold(X=X, y=y_shuff, PMA_scan=PMA_scan, 
                               corr_idx=for_corr_index, n_folds=n_folds, mode='cont')
            rho, pval, mae, r2score = evaluate(y_true=y_true, y_pred=y_pred, mode='cont')
            res['rho'].append(rho)
            res['pval'].append(pval)
            res['mae'].append(mae)
            res['r2'].append(r2score)		
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
        res = {'rho' : [], 'pval' : [], 'mae' : [], 'r2' : []}
        for i in range(num_reps):
            X,y,PMA_scan = unison_shuffled_copies(a=X, b=y, c=PMA_scan)
            y_true, y_pred = run_Kfold(X=X, y=y, PMA_scan=PMA_scan, 
                               corr_idx=for_corr_index, n_folds=n_folds, mode='cont')
            rho, pval, mae, r2score = evaluate(y_true=y_true, y_pred=y_pred, mode='cont')
            res['rho'].append(rho)
            res['pval'].append(pval)
            res['mae'].append(mae)
            res['r2'].append(r2score)
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
                           mode=mode)
        for el in range(len(y_test)):
                y_true.append(y_test[el])
                y_pred.append(y_out[el]) 
    return y_true, y_pred
        
def evaluate(y_true, y_pred, mode):
    if mode == 'cat':
        ## auc, balanced acc, specificity, sensitivity 
        return roc_auc_score(y_true, y_pred), balanced_accuracy_score(y_true, y_pred), recall_score(y_true, y_pred, pos_label=0),recall_score(y_true, y_pred, pos_label=1) 
    if mode == 'cont':
        #print('True:', y_true)
        #print('Predicted:', y_pred)
        return pg.corr(y_true, y_pred)['r'][0], pg.corr(y_true, y_pred)['p-val'][0], mean_absolute_error(y_true, y_pred), r2_score(y_true, y_pred)
        

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
o_file = '../../Results/predictions/minAD_outcomes_prediction_results.txt'
if os.path.exists(o_file):
    os.remove(o_file)

### stupidest ever but I cannot survive another loop 
# ===================== COGNITIVE =======#
outcome='Cognitive Score'
for cohort in ['A', 'C', 'D', 'E']:

	with open(o_file, 'a') as the_file:
		the_file.write('Cohort {}\n'.format(cohort))
		the_file.write('{}\n'.format(outcome))
		the_file.write('Inputs, RHO mean (std), Pval mean (std), MAE mean (std), Permuted RHO (p-val), Permuted MAE (p-val)\n')

	if cohort == 'E':		
		sub_df = df[df.subject_id.isin(cohorts['D'])].copy()
		sub_df = sub_df.sample(frac=1).reset_index(drop=True)
		sub_df = sub_df.head(n_folds).copy()	

	else: 
		sub_df = df[df.subject_id.isin(cohorts[cohort])].copy()
		sub_df = sub_df.sample(frac=1).reset_index(drop=True)
	
	print('Cohort {}'.format(cohort))
	print('{}'.format(outcome))
	print('Inputs, RHO mean (std), Pval mean (std), MAE mean (std), Permuted RHO (p-val), Permuted MAE (p-val)')

	for key in inputs.keys():
    
		in_cols, for_corr_index = get_columns_to_correct(inDict=inputs, key=key)
    
    		### get inputs
		X = sub_df[in_cols].values
		y = sub_df[outcome].values
		PMA_scan = sub_df['PMA_scan'].values
    
    		### RUN & EVALUATE 
		res = get_repeated_Kfold(X=X, y=y, PMA_scan=PMA_scan, 
                             corr_idx=for_corr_index, n_folds=n_folds, 
                             mode='cont', num_reps=num_rep)
		null = get_null_distribution(X=X, y=y, PMA_scan=PMA_scan, 
                             corr_idx=for_corr_index, n_folds=n_folds, 
                             mode='cont', it=it_num)

		rho_p = np.sum(null['rho']>=res['rho'][0]) / len(null['rho'])
		mae_p = np.sum(null['mae']<=res['mae'][0]) / len(null['mae'])

		out_message = '{}, {:.3f}({:.3f}),{:.3f}({:.3f}),{:.3f}({:.3f}),{:.3f}({:.3f}),{:.3f}({:.6f}),{:.3f}({:.6f})\n'.format(
                        key,
                        np.mean(res['rho']), np.std(res['rho']),
                        np.mean(res['pval']), np.std(res['pval']),
                        np.mean(res['mae']), np.std(res['mae']),
			 np.mean(res['r2']), np.std(res['r2']),
                        res['rho'][0],rho_p,
                        res['mae'][0],mae_p
    			)
		print(out_message)
		with open(o_file, 'a') as the_file:
    			the_file.write(out_message)
	
	for random_parc_num in [128, 256, 512, 1024]:
		path_to_random = '/neurospin/grip/external_databases/dHCP_CR_JD_2018/Projects/eLife_replication/DerivedData/extracted_metrics/'
		random_df = pd.read_csv(os.path.join(path_to_random, 'random_parcellation_{}_FA_median_minAD.csv'.format(random_parc_num)), index_col=0)
		random_df = pd.merge(random_df, sub_df[['subject_id', 'PMA_scan','Cognitive Score','Language Score','Motor Score']], on='subject_id', how='inner')

		random_FA = [col for col in random_df.columns if 'FA' in col]
		random_inputs = {
				 #'Random {}'.format(random_parc_num) : [random_FA, np.zeros_like(random_FA, dtype=int)], 
				 'Random {} corrected'.format(random_parc_num) : [random_FA, np.ones_like(random_FA, dtype=int)], 
				} 

		if cohort == 'E':
			random_df = random_df[random_df.subject_id.isin(cohorts['D'])].copy()
			## shuffle the data 
			random_df = random_df.sample(frac=1).reset_index(drop=True)
			random_df = random_df.head(n_folds).copy()
		else: 
			random_df = random_df[random_df.subject_id.isin(cohorts[cohort])].copy()
			## shuffle the data 
			random_df = random_df.sample(frac=1).reset_index(drop=True)

		for key in random_inputs.keys():
    
			in_cols, for_corr_index = get_columns_to_correct(inDict=random_inputs, key=key)
    
    			### get inputs
			X = random_df[in_cols].values
			y = random_df[outcome].values
			PMA_scan = random_df['PMA_scan'].values
    
    			### RUN & EVALUATE 
			res = get_repeated_Kfold(X=X, y=y, PMA_scan=PMA_scan, 
                             corr_idx=for_corr_index, n_folds=n_folds, 
                             mode='cont', num_reps=num_rep)
			null = get_null_distribution(X=X, y=y, PMA_scan=PMA_scan, 
                             corr_idx=for_corr_index, n_folds=n_folds, 
                             mode='cont', it=it_num)
    
			rho_p = np.sum(null['rho']>=res['rho'][0]) / len(null['rho'])
			mae_p = np.sum(null['mae']<=res['mae'][0]) / len(null['mae'])

			out_message = '{}, {:.3f}({:.3f}),{:.3f}({:.3f}),{:.3f}({:.3f}),{:.3f}({:.3f}),{:.3f}({:.6f}),{:.3f}({:.6f})\n'.format(
                        		key,
                        	np.mean(res['rho']), np.std(res['rho']),
                        	np.mean(res['pval']), np.std(res['pval']),
                        	np.mean(res['mae']), np.std(res['mae']),
				 np.mean(res['r2']), np.std(res['r2']),
                        	res['rho'][0],rho_p,
                        	res['mae'][0],mae_p
    				)
			print(out_message)
			with open(o_file, 'a') as the_file:
    				the_file.write(out_message)



# ===================== Language =======#
outcome='Language Score'
for cohort in ['A', 'C', 'D', 'E']:

	with open(o_file, 'a') as the_file:
		the_file.write('Cohort {}\n'.format(cohort))
		the_file.write('{}\n'.format(outcome))
		the_file.write('Inputs, RHO mean (std), Pval mean (std), MAE mean (std), Permuted RHO (p-val), Permuted MAE (p-val)\n')

	if cohort == 'E':		
		sub_df = df[df.subject_id.isin(cohorts['D'])].copy()
		sub_df = sub_df.sample(frac=1).reset_index(drop=True)
		sub_df = sub_df.head(n_folds).copy()	

	else: 
		sub_df = df[df.subject_id.isin(cohorts[cohort])].copy()
		sub_df = sub_df.sample(frac=1).reset_index(drop=True)
	
	print('Cohort {}'.format(cohort))
	print('{}'.format(outcome))
	print('Inputs, RHO mean (std), Pval mean (std), MAE mean (std), Permuted RHO (p-val), Permuted MAE (p-val)')

	for key in inputs.keys():
    
		in_cols, for_corr_index = get_columns_to_correct(inDict=inputs, key=key)
    
    		### get inputs
		X = sub_df[in_cols].values
		y = sub_df[outcome].values
		PMA_scan = sub_df['PMA_scan'].values
    
    		### RUN & EVALUATE 
		res = get_repeated_Kfold(X=X, y=y, PMA_scan=PMA_scan, 
                             corr_idx=for_corr_index, n_folds=n_folds, 
                             mode='cont', num_reps=num_rep)
		null = get_null_distribution(X=X, y=y, PMA_scan=PMA_scan, 
                             corr_idx=for_corr_index, n_folds=n_folds, 
                             mode='cont', it=it_num)

		rho_p = np.sum(null['rho']>=res['rho'][0]) / len(null['rho'])
		mae_p = np.sum(null['mae']<=res['mae'][0]) / len(null['mae'])

		out_message = '{}, {:.3f}({:.3f}),{:.3f}({:.3f}),{:.3f}({:.3f}),{:.3f}({:.3f}),{:.3f}({:.6f}),{:.3f}({:.6f})\n'.format(
                        key,
                        np.mean(res['rho']), np.std(res['rho']),
                        np.mean(res['pval']), np.std(res['pval']),
                        np.mean(res['mae']), np.std(res['mae']),
			 np.mean(res['r2']), np.std(res['r2']),
                        res['rho'][0],rho_p,
                        res['mae'][0],mae_p
    			)
		print(out_message)
		with open(o_file, 'a') as the_file:
    			the_file.write(out_message)
	
	for random_parc_num in [128, 256, 512, 1024]:
		path_to_random = '/neurospin/grip/external_databases/dHCP_CR_JD_2018/Projects/eLife_replication/DerivedData/extracted_metrics/'
		random_df = pd.read_csv(os.path.join(path_to_random, 'random_parcellation_{}_FA_median_minAD.csv'.format(random_parc_num)), index_col=0)
		random_df = pd.merge(random_df, sub_df[['subject_id', 'PMA_scan','Cognitive Score','Language Score','Motor Score']], on='subject_id', how='inner')

		random_FA = [col for col in random_df.columns if 'FA' in col]
		random_inputs = {
				 #'Random {}'.format(random_parc_num) : [random_FA, np.zeros_like(random_FA, dtype=int)], 
				 'Random {} corrected'.format(random_parc_num) : [random_FA, np.ones_like(random_FA, dtype=int)], 
				} 

		if cohort == 'E':
			random_df = random_df[random_df.subject_id.isin(cohorts['D'])].copy()
			## shuffle the data 
			random_df = random_df.sample(frac=1).reset_index(drop=True)
			random_df = random_df.head(n_folds).copy()
		else: 
			random_df = random_df[random_df.subject_id.isin(cohorts[cohort])].copy()
			## shuffle the data 
			random_df = random_df.sample(frac=1).reset_index(drop=True)

		for key in random_inputs.keys():
    
			in_cols, for_corr_index = get_columns_to_correct(inDict=random_inputs, key=key)
    
    			### get inputs
			X = random_df[in_cols].values
			y = random_df[outcome].values
			PMA_scan = random_df['PMA_scan'].values
    
    			### RUN & EVALUATE 
			res = get_repeated_Kfold(X=X, y=y, PMA_scan=PMA_scan, 
                             corr_idx=for_corr_index, n_folds=n_folds, 
                             mode='cont', num_reps=num_rep)
			null = get_null_distribution(X=X, y=y, PMA_scan=PMA_scan, 
                             corr_idx=for_corr_index, n_folds=n_folds, 
                             mode='cont', it=it_num)
    
			rho_p = np.sum(null['rho']>=res['rho'][0]) / len(null['rho'])
			mae_p = np.sum(null['mae']<=res['mae'][0]) / len(null['mae'])

			out_message = '{}, {:.3f}({:.3f}),{:.3f}({:.3f}),{:.3f}({:.3f}),{:.3f}({:.3f}),{:.3f}({:.6f}),{:.3f}({:.6f})\n'.format(
                        		key,
                        	np.mean(res['rho']), np.std(res['rho']),
                        	np.mean(res['pval']), np.std(res['pval']),
                        	np.mean(res['mae']), np.std(res['mae']),
				 np.mean(res['r2']), np.std(res['r2']),
                        	res['rho'][0],rho_p,
                        	res['mae'][0],mae_p
    				)
			print(out_message)
			with open(o_file, 'a') as the_file:
    				the_file.write(out_message)

# ===================== Motor =======#
outcome='Motor Score'
for cohort in ['A', 'C', 'D', 'E']:

	with open(o_file, 'a') as the_file:
		the_file.write('Cohort {}\n'.format(cohort))
		the_file.write('{}\n'.format(outcome))
		the_file.write('Inputs, RHO mean (std), Pval mean (std), MAE mean (std), Permuted RHO (p-val), Permuted MAE (p-val)\n')

	if cohort == 'E':		
		sub_df = df[df.subject_id.isin(cohorts['D'])].copy()
		sub_df = sub_df.sample(frac=1).reset_index(drop=True)
		sub_df = sub_df.head(n_folds).copy()	

	else: 
		sub_df = df[df.subject_id.isin(cohorts[cohort])].copy()
		sub_df = sub_df.sample(frac=1).reset_index(drop=True)
	
	print('Cohort {}'.format(cohort))
	print('{}'.format(outcome))
	print('Inputs, RHO mean (std), Pval mean (std), MAE mean (std), Permuted RHO (p-val), Permuted MAE (p-val)')

	for key in inputs.keys():
    
		in_cols, for_corr_index = get_columns_to_correct(inDict=inputs, key=key)
    
    		### get inputs
		X = sub_df[in_cols].values
		y = sub_df[outcome].values
		PMA_scan = sub_df['PMA_scan'].values
    
    		### RUN & EVALUATE 
		res = get_repeated_Kfold(X=X, y=y, PMA_scan=PMA_scan, 
                             corr_idx=for_corr_index, n_folds=n_folds, 
                             mode='cont', num_reps=num_rep)
		null = get_null_distribution(X=X, y=y, PMA_scan=PMA_scan, 
                             corr_idx=for_corr_index, n_folds=n_folds, 
                             mode='cont', it=it_num)

		rho_p = np.sum(null['rho']>=res['rho'][0]) / len(null['rho'])
		mae_p = np.sum(null['mae']<=res['mae'][0]) / len(null['mae'])

		out_message = '{}, {:.3f}({:.3f}),{:.3f}({:.3f}),{:.3f}({:.3f}),{:.3f}({:.3f}),{:.3f}({:.6f}),{:.3f}({:.6f})\n'.format(
                        key,
                        np.mean(res['rho']), np.std(res['rho']),
                        np.mean(res['pval']), np.std(res['pval']),
                        np.mean(res['mae']), np.std(res['mae']),
			 np.mean(res['r2']), np.std(res['r2']),
                        res['rho'][0],rho_p,
                        res['mae'][0],mae_p
    			)
		print(out_message)
		with open(o_file, 'a') as the_file:
    			the_file.write(out_message)
	
	for random_parc_num in [128, 256, 512, 1024]:
		path_to_random = '/neurospin/grip/external_databases/dHCP_CR_JD_2018/Projects/eLife_replication/DerivedData/extracted_metrics/'
		random_df = pd.read_csv(os.path.join(path_to_random, 'random_parcellation_{}_FA_median_minAD.csv'.format(random_parc_num)), index_col=0)
		random_df = pd.merge(random_df, sub_df[['subject_id', 'PMA_scan','Cognitive Score','Language Score','Motor Score']], on='subject_id', how='inner')

		random_FA = [col for col in random_df.columns if 'FA' in col]
		random_inputs = {
				 #'Random {}'.format(random_parc_num) : [random_FA, np.zeros_like(random_FA, dtype=int)], 
				 'Random {} corrected'.format(random_parc_num) : [random_FA, np.ones_like(random_FA, dtype=int)], 
				} 

		if cohort == 'E':
			random_df = random_df[random_df.subject_id.isin(cohorts['D'])].copy()
			## shuffle the data 
			random_df = random_df.sample(frac=1).reset_index(drop=True)
			random_df = random_df.head(n_folds).copy()
		else: 
			random_df = random_df[random_df.subject_id.isin(cohorts[cohort])].copy()
			## shuffle the data 
			random_df = random_df.sample(frac=1).reset_index(drop=True)

		for key in random_inputs.keys():
    
			in_cols, for_corr_index = get_columns_to_correct(inDict=random_inputs, key=key)
    
    			### get inputs
			X = random_df[in_cols].values
			y = random_df[outcome].values
			PMA_scan = random_df['PMA_scan'].values
    
    			### RUN & EVALUATE 
			res = get_repeated_Kfold(X=X, y=y, PMA_scan=PMA_scan, 
                             corr_idx=for_corr_index, n_folds=n_folds, 
                             mode='cont', num_reps=num_rep)
			null = get_null_distribution(X=X, y=y, PMA_scan=PMA_scan, 
                             corr_idx=for_corr_index, n_folds=n_folds, 
                             mode='cont', it=it_num)
    
			rho_p = np.sum(null['rho']>=res['rho'][0]) / len(null['rho'])
			mae_p = np.sum(null['mae']<=res['mae'][0]) / len(null['mae'])

			out_message = '{}, {:.3f}({:.3f}),{:.3f}({:.3f}),{:.3f}({:.3f}),{:.3f}({:.3f}),{:.3f}({:.6f}),{:.3f}({:.6f})\n'.format(
                        		key,
                        	np.mean(res['rho']), np.std(res['rho']),
                        	np.mean(res['pval']), np.std(res['pval']),
                        	np.mean(res['mae']), np.std(res['mae']),
				 np.mean(res['r2']), np.std(res['r2']),
                        	res['rho'][0],rho_p,
                        	res['mae'][0],mae_p
    				)
			print(out_message)
			with open(o_file, 'a') as the_file:
    				the_file.write(out_message)

