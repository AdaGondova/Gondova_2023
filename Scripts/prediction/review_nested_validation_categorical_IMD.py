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
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, balanced_accuracy_score, mean_absolute_error, recall_score
import statsmodels.api as sm

#### ========================= GLOBAL ============================= ###
with open(r"../../DerivedData/cohorts_subjects_list.pickle", "rb") as input_file:
        cohorts = pickle.load(input_file)

outcomes = ['group_cat','Cognitive_cat','Language_cat','Motor_cat']
n_folds = 46 # to recreate independent validation set
opt_itr = 5 # inner loop 
N = 5 # outer loop 

param_grid = {
    
    #'nu' : np.linspace(0.1,1,10),
    'C' : [0.01,0.1,1,5,9,10],
    'kernel' : ['linear'],
    'shrinking' : [True, False]
             }
#### ======================== FUNCTIONS ============== #
def evaluate(y_true, y_pred):
    	     
    return roc_auc_score(y_true, y_pred), balanced_accuracy_score(y_true, y_pred), recall_score(y_true, y_pred, pos_label=0),recall_score(y_true, y_pred, pos_label=1) 


def preprocess(X_train, X_val, cols, inflection=36, correct_motion=False):
    ## inpute
    X_train, X_val = _inpute_median(X_train=X_train, X_val=X_val, cols=train_cols)
    
    #if correct_motion == True:
    #    X_train, X_val = _correct_motion(X_train=X_train, X_val=X_val, cols=train_cols)
    ## correct PMA
    #X_train, X_val = _correct_PMA_scan(X_train=X_train, X_val=X_val, cols=train_cols, inflection=36)
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
df = pd.read_csv('../for_risk_and_env_factors.csv', index_col=0)


df['group'] = 'FT'
df.loc[df['GA_birth'] < 37, 'group'] = 'PT'

df['group_cat'] = 0
df.loc[df['GA_birth'] < 37, 'group_cat'] = 1

scores = ['Cognitive','Language','Motor']
for score in scores: 
	df[score+'_cat'] = 0
	df.loc[df[score] <= 85, score+'_cat'] = 1


### ======================= TRAIN/EVALUATE ============================ ###
## hyperparameters
param_list = list(ParameterSampler(param_grid, n_iter=opt_itr, random_state=42))
#for d in param_list:
#    d['nu'] = np.round(d['nu'], 1)

outer_res = []
df = df[df.subject_id.isin(cohorts['A'])].copy()
for inputs in ['IMD', 'IMD_and_age']:

	if inputs == 'IMD':
		train_cols = ['corr_FA', 'sex', 'IMD']
	else:
		train_cols = ['corr_FA', 'GA_birth', 'sex', 'IMD']
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
					model = SVC(**rep)
					model.fit(X_train[train_cols].values, X_train[outcome].values)
					y_out = model.predict(X_val[train_cols].values)
        
					y_hat.append(y_out[0])
					y_true.append(X_val[outcome].values[0])
       
				AUC, ACC, SPEC, SENS = evaluate(y_true=y_true, y_pred=y_hat)
				inner_res.append([i, AUC, ACC, SPEC, SENS])
				print('{}th Inner Loop FINISHED'.format(i))
				print(AUC, ACC, SPEC, SENS)
        	
        
			inner_res = pd.DataFrame(data=inner_res, columns=['params', 'AUC', 'ACC', 'SPEC', 'SENS'])
			opt_params = param_list[inner_res[inner_res['AUC'] == inner_res['AUC'].max()]['params'].values[0]]   

			## outer loop 
			## preprocess 
			df_train, df_test = preprocess(X_train=df_train, X_val=df_test, cols=train_cols, inflection=36, correct_motion=True)

			fmodel = SVC(**opt_params)
			fmodel.fit(df_train[train_cols].values, df_train[outcome].values)
			f_yhat = fmodel.predict(df_test[train_cols].values).round(0)
			AUC, ACC, SPEC, SENS = evaluate(y_true=df_test[outcome].values, y_pred=f_yhat)
			
			
			outer_res.append([outcome,inputs, opt_params,AUC, ACC, SPEC, SENS, np.array(df_test[outcome].values),np.array(f_yhat)])
    
			print('\n***')
			print('{}th Outer Loop FINISHED'.format(n))
			print(AUC, ACC, SPEC, SENS)
			print('\n***')
    
outer_res = pd.DataFrame(data=outer_res, columns=['outcome','inputs', 'params', 'AUC', 'ACC', 'SPEC', 'SENS', 'y_true', 'y_pred'])  
outer_res.to_csv('../../Results/predictions/review_validation_categortical_IMD.csv')



	

