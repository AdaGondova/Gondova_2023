OUTDATED_IGNORE=1

import pandas as pd 
import numpy as np 
from sklearn.svm import NuSVR, NuSVC
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import LeaveOneOut, KFold, StratifiedKFold, train_test_split, GridSearchCV, ParameterSampler
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
import statsmodels.api as sm

#### ========================= GLOBAL ============================= ###
with open(r"../../DerivedData/cohorts_subjects_list.pickle", "rb") as input_file:
        cohorts = pickle.load(input_file)

#outcomes = ['Cognitive Score','Language Score','Motor Score']
outcomes = ['GA_birth']
n_folds = 46 # to recreate independent validation set
opt_itr = 4 # inner loop 
N = 5 # outer loop 

param_grid = {
    
    'nu' : np.linspace(0.1,1,10),
    'C' : [0.01,0.1,1,5,9,10],
    'kernel' : ['linear'],
    'shrinking' : [True, False]
             }
#### ======================== FUNCTIONS ============== #
def evaluate(y_true, y_pred):
    rho = pg.corr(y_true, y_pred)['r'][0]
    if np.isnan(rho):
        rho = 'undefined'	     
    return rho, mean_absolute_error(y_true, y_pred), r2_score(y_true, y_pred)


def preprocess(X_train, X_val, cols, inflection=36, correct_motion=False):
    ## inpute
    X_train, X_val = _inpute_median(X_train=X_train, X_val=X_val, cols=train_cols)
    
    if correct_motion == True:
        X_train, X_val = _correct_motion(X_train=X_train, X_val=X_val, cols=train_cols)
    ## correct PMA
    X_train, X_val = _correct_PMA_scan(X_train=X_train, X_val=X_val, cols=train_cols, inflection=36)
    ## scale
    X_train, X_val = _scaling(X_train=X_train, X_val=X_val, cols=train_cols)
    
    return X_train, X_val

def _inpute_median(X_train, X_val, cols):
    for col in cols:
        md = np.nanmedian(X_train[col])
        
        X_train[col].fillna(md, inplace= True)
        X_val[col].fillna(md, inplace= True)
   
    return X_train, X_val

def _correct_PMA_scan(X_train, X_val, cols, inflection = 36):
    
    x0 = np.array([min(X_train.PMA_scan.values), inflection, max(X_train.PMA_scan.values)])
    
    for col in cols:
        
        myPWLF = pwlf.PiecewiseLinFit(X_train.PMA_scan.values, X_train[col].values)
        myPWLF.fit_with_breaks(x0)
        
        ## correct train 
        yHat_train = myPWLF.predict(X_train.PMA_scan.values)
        res_train = X_train[col].values - yHat_train
        
        ## correct test 
        yHat_test = myPWLF.predict(X_val.PMA_scan.values)
        res_test = X_val[col].values - yHat_test
        
        X_train[col] = res_train
        X_val[col] = res_test
    return X_train, X_val   

def _scaling(X_train, X_val, cols):
    
    scaler = MinMaxScaler()
    scaler.fit(X_train[cols].values)
    
    X_train[cols] = scaler.transform(X_train[cols].values)
    X_val[cols] = scaler.transform(X_val[cols].values)
    
    return X_train, X_val

def _correct_motion(X_train, X_val, cols):
    
    for col in cols:
        
        model = sm.OLS(X_train[col].values, X_train[['qc_translation', 'qc_rotation']].values).fit()
        
        yHat_train = model.predict(X_train[['qc_translation', 'qc_rotation']].values) 
        res_train = X_train[col].values - yHat_train
        #print(yHat_train)
    
    
        yHat_test = model.predict(X_val[['qc_translation', 'qc_rotation']].values) 
        res_test = X_val[col].values - yHat_test
        
        X_train[col] = res_train
        X_val[col] = res_test
    return X_train, X_val   
#### ========================= READ IN DATA ============================= ###### ages
df = pd.read_csv('../../DerivedData/cohortA_subjects_clinical.csv', index_col=0)

### global FA
glob = pd.read_csv('../../DerivedData/extracted_metrics/global_cortical_diffusion_metrics_median.csv', index_col=0)

metrics = ['FA', 'L1', 'RD', 'MD']
hemispheres = ['left', 'right']
for metric in metrics:
    for i, row in glob.iterrows():
        glob.loc[i,metric] = np.mean([row['left_{}'.format(metric)], 
                                   row['right_{}'.format(metric)]])
df = pd.merge(df, glob[['subject_id', 'FA']], on=['subject_id'])  

### median FA regions 
diff = pd.read_csv('../../DerivedData/extracted_metrics/neonat_segmentation_diffusion_metric_median.csv', index_col=0)

FA_cols = [col for col in diff.columns if 'FA' in col]
FA_cols.extend(['subject_id'])

#new_df = df.copy()
df = pd.merge(df, diff[FA_cols], on=['subject_id'])
df = df[df.subject_id.isin(cohorts['D'])]

### ADD MOTION VALUES HERE FOR THE CORRECTION!!! 
qc = pd.read_csv('../../SourceData/release3_subject_info.tsv', sep='\t')
for i, row in df.iterrows():
    
    trans = qc[(qc['participant_id '] == row.subject_id + ' ') & (qc['session_id '] == row.session_id)]['qc_dmri_shard_translation '].values[0]
    rot = qc[(qc['participant_id '] == row.subject_id + ' ') & (qc['session_id '] == row.session_id)]['qc_dmri_shard_rotation '].values[0]
    #snr = qc[(qc['participant_id '] == row.subject_id + ' ') & (qc['session_id '] == row.session_id)]['qc_dmri_shard_snr '].values[0]
    #outlier = qc[(qc['participant_id '] == row.subject_id + ' ') & (qc['session_id '] == row.session_id)]['qc_dmri_shard_outlier_ratio '].values[0]
    
    df.loc[i, 'qc_translation'] = trans
    df.loc[i, 'qc_rotation'] = rot
    #df.loc[i, 'qc_dmri_shard_snr'] = snr
    #df.loc[i, 'qc_dmri_shard_outlier_ratio'] = outlier

### ======================= TRAIN/EVALUATE ============================ ###
## hyperparameters
param_list = list(ParameterSampler(param_grid, n_iter=opt_itr, random_state=42))
for d in param_list:
    d['nu'] = np.round(d['nu'], 1)

outer_res = []

for inputs in ['ROIs age corrected', 128, 256, 512, 1024]:
	print(inputs)
	if inputs == 'ROIs age corrected':
		train_cols = FA_cols[:-1]


	else: 
		path_to_random = '/neurospin/grip/external_databases/dHCP_CR_JD_2018/Projects/eLife_replication/DerivedData/extracted_metrics/'
		random_df = pd.read_csv(os.path.join(path_to_random, 'random_parcellation_{}_diffusion_metric_median.csv'.format(inputs)), index_col=0)
		df = pd.merge(random_df, df[['subject_id', 'PMA_scan','Cognitive Score','Language Score','Motor Score','GA_birth', 'qc_translation', 'qc_rotation']], on='subject_id', how='inner').copy()
		train_cols = [col for col in random_df.columns if 'FA' in col]


	### loop over outcomes:
	for outcome in outcomes:
		print(outcome)
		for n in range(N):
    
			df_train, df_test = train_test_split( df,  test_size=n_folds, random_state=n)
			df_train.reset_index(inplace=True, drop=True)
			df_test.reset_index(inplace=True, drop=True)

			### inner loop
			inner_res = []
			loo = LeaveOneOut()
			for i,  rep in enumerate(param_list):    
				y_hat = []
				y_true = []

				for train_index, test_index in loo.split(df_train):
        
					X_train, X_val = df_train.loc[train_index], df_train.loc[test_index]
					#preprocess 
					X_train, X_val = preprocess(X_train=X_train, X_val=X_val, cols=train_cols, inflection=36, correct_motion=True)
    
					#
					model = NuSVR(**rep)
					model.fit(X_train[train_cols].values, X_train[outcome].values)
					y_out = model.predict(X_val[train_cols].values)
        
					y_hat.append(y_out[0])
					y_true.append(X_val[outcome].values[0])
        
				rho, mae, r2 = evaluate(y_true=y_true, y_pred=y_hat)
				inner_res.append([i, rho, mae, r2])
				print('{}th Inner Loop FINISHED'.format(i))
				print(rho, mae, r2)
        	
        
			inner_res = pd.DataFrame(data=inner_res, columns=['params', 'rho', 'mae', 'r2'])
			opt_params = param_list[inner_res[inner_res['mae'] == inner_res['mae'].min()]['params'].values[0]]   

			## outer loop 
			## preprocess 
			df_train, df_test = preprocess(X_train=df_train, X_val=df_test, cols=train_cols, inflection=36, correct_motion=True)

			fmodel = NuSVR(**opt_params)
			fmodel.fit(df_train[train_cols].values, df_train[outcome].values)
			f_yhat = fmodel.predict(df_test[train_cols].values).round(0)
			rho, mae, r2 = evaluate(y_true=df_test[outcome].values, y_pred=f_yhat)
			
			
			outer_res.append([outcome,inputs, opt_params,rho, mae, r2, np.array(df_test[outcome].values),np.array(f_yhat)])
    
			print('\n***')
			print('{}th Outer Loop FINISHED'.format(n))
			print(rho, mae, r2)
			print('\n***')
    
outer_res = pd.DataFrame(data=outer_res, columns=['outcome','inputs', 'params', 'rho', 'mae', 'r2', 'y_true', 'y_pred'])  
outer_res.to_csv('../../Results/predictions/nested_validation_GA_birth_motion_corrected.csv')

outer_res = []

for inputs in ['ROIs age corrected', 128, 256, 512, 1024]:
	print(inputs)
	if inputs == 'ROIs age corrected':
		train_cols = FA_cols[:-1]


	else: 
		path_to_random = '/neurospin/grip/external_databases/dHCP_CR_JD_2018/Projects/eLife_replication/DerivedData/extracted_metrics/'
		random_df = pd.read_csv(os.path.join(path_to_random, 'random_parcellation_{}_diffusion_metric_median.csv'.format(inputs)), index_col=0)
		df = pd.merge(random_df, df[['subject_id', 'PMA_scan','Cognitive Score','Language Score','Motor Score','GA_birth', 'qc_translation', 'qc_rotation']], on='subject_id', how='inner').copy()
		train_cols = [col for col in random_df.columns if 'FA' in col]


	### loop over outcomes:
	for outcome in outcomes:
		print(outcome)
		for n in range(N):
    
			df_train, df_test = train_test_split( df,  test_size=n_folds, random_state=n)
			df_train.reset_index(inplace=True, drop=True)
			df_test.reset_index(inplace=True, drop=True)

			### inner loop
			inner_res = []
			loo = LeaveOneOut()
			for i,  rep in enumerate(param_list):    
				y_hat = []
				y_true = []

				for train_index, test_index in loo.split(df_train):
        
					X_train, X_val = df_train.loc[train_index], df_train.loc[test_index]
					#preprocess 
					X_train, X_val = preprocess(X_train=X_train, X_val=X_val, cols=train_cols, inflection=36, correct_motion=False)
    
					#
					model = NuSVR(**rep)
					model.fit(X_train[train_cols].values, X_train[outcome].values)
					y_out = model.predict(X_val[train_cols].values)
        
					y_hat.append(y_out[0])
					y_true.append(X_val[outcome].values[0])
        
				rho, mae, r2 = evaluate(y_true=y_true, y_pred=y_hat)
				inner_res.append([i, rho, mae, r2])
				print('{}th Inner Loop FINISHED'.format(i))
				print(rho, mae, r2)
        	
        
			inner_res = pd.DataFrame(data=inner_res, columns=['params', 'rho', 'mae', 'r2'])
			opt_params = param_list[inner_res[inner_res['mae'] == inner_res['mae'].min()]['params'].values[0]]   

			## outer loop 
			## preprocess 
			df_train, df_test = preprocess(X_train=df_train, X_val=df_test, cols=train_cols, inflection=36, correct_motion=False)

			fmodel = NuSVR(**opt_params)
			fmodel.fit(df_train[train_cols].values, df_train[outcome].values)
			f_yhat = fmodel.predict(df_test[train_cols].values).round(0)
			rho, mae, r2 = evaluate(y_true=df_test[outcome].values, y_pred=f_yhat)
			
			
			outer_res.append([outcome,inputs, opt_params,rho, mae, r2, np.array(df_test[outcome].values),np.array(f_yhat)])
    
			print('\n***')
			print('{}th Outer Loop FINISHED'.format(n))
			print(rho, mae, r2)
			print('\n***')
    
outer_res = pd.DataFrame(data=outer_res, columns=['outcome','inputs', 'params', 'rho', 'mae', 'r2', 'y_true', 'y_pred'])  
outer_res.to_csv('../../Results/predictions/nested_validation_GA_birth_age_corrected.csv')

	

