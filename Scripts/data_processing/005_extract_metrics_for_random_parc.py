from soma import aims 
import matplotlib.pyplot as plt 
import numpy as np 
import pandas as pd
import os, datetime, sys

### requires bv 
### I am adjusting this slightly to get the median as well - might be more reliable regional descriptor

def read_texture(path):
    text_f = aims.read(path)
    text_a = np.array(text_f)[0]
    
    return text_a



if __name__ == "__main__":
	print('START at: {}'.format(datetime.datetime.now()))

	### read in the file containing the list of subject ID and session IDs
	if len(sys.argv) < 2:
		print("You must provide subject file!")
		sys.exit()
	else:
		subject_file = sys.argv[1]

	subjects = pd.read_csv(subject_file, names=['subject_id', 'session_id', 'template'])

	
	for segmentation in [64, 128, 256, 512]:
		print('Working on {} extraction'.format(segmentation))
    		df = pd.DataFrame()
    		for i, row in subjects.iterrows():   
        		if i % 100 == 0 :
            			print('Running {}/{}'.format(i,len(subjects)))
        		
			for hemi in ['left', 'right']:
            
            			iDIR = '../../DerivedData/subjects/sub-{}/ses-{}'.format(row.subject_id, row.session_id)        
            			iFile = 'sub-{}_ses-{}_hemi-{}_pathKmeans_{}.label.gii'.format(row.subject_id, row.session_id, hemi, segmentation )
        
            			iSeg = os.path.join(iDIR, iFile)
              
            			if os.path.isfile(iSeg):
                
                			segm = read_texture(iSeg)
                			labels = np.unique(segm)[:-1]  ## largest label == cingulate - remove 

                			for metric in ['FA', 'L1', 'RD', 'MD']:
                    
                    				iMetric = os.path.join(iDIR, 'sub-{}_ses-{}_{}_texture_{}_majority.gii'.format(row.subject_id, row.session_id, 
                                                                                                   hemi, metric))
                    				if os.path.isfile(iMetric):
                      	 				met = read_texture(iMetric)
                
                        				df.loc[i, 'subject_id'] = row.subject_id
                        				df.loc[i, 'session_id'] = row.session_id
                
                					
                        				for label in labels:
                            					#df.loc[i, '{}_{}_{}'.format(hemi, label, metric)] = np.mean(met[segm==label].ravel()) 
								df.loc[i, '{}_{}_{}'.format(hemi, label, metric)] = np.median(met[segm==label].ravel()) 
                   	 			else: 
                       		 			print(row.subject_id, row.session_id, '{} file missing'.format(metric))
            			else:
                			print(row.subject_id, row.session_id, 'segmentation missing')
                
    		df['session_id'] = df['session_id'].astype(int)
    		#df.to_csv('../../DerivedData/extracted_metrics/random_parcellation_{}_diffusion_metric.csv'.format(segmentation * 2))
		df.to_csv('../../DerivedData/extracted_metrics/random_parcellation_{}_diffusion_metric_median.csv'.format(segmentation * 2))

		#if os.path.isfile('../../DerivedData/extracted_metrics/random_parcellation_{}_diffusion_metric.csv'):
		if os.path.isfile('../../DerivedData/extracted_metrics/random_parcellation_{}_diffusion_metric_median.csv'):
			print('{} finished'.format(segmentation*2))
	
	print('END at: {}'.format(datetime.datetime.now()))
