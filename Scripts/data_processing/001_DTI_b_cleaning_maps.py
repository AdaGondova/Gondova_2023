## requires bv bash for the soma.aims functionality

from soma import aims 
import numpy as np 
import os, sys, datetime
import pandas as pd

# ========================= #

def correct_DTI(subj_id, ses_id):

    incorrectIdx = []

    for metric in ['FA', 'L1', 'L2', 'L3']:

        vol = aims.read('/neurospin/grip/external_databases/dHCP_CR_JD_2018/release3/dhcp_dmri_shard_pipeline/sub-{}/ses-{}/dwi/DTI/dtifit_b1000/sub-{}_ses-{}_{}.nii.gz'.format(
                                    subj_id, ses_id,subj_id, ses_id, metric))
        im = vol.arraydata()[0]
    
        args = np.argwhere(im < 0)
        for arg in args:    
            incorrectIdx.append(arg)
        
        if metric == 'FA':
            args = np.argwhere(im > 1)
        for arg in args:    
            incorrectIdx.append(arg)
        
    if len(incorrectIdx) != 0:
        new = [tuple(row) for row in incorrectIdx]
        uniques = np.unique(new, axis=0)


        for metric in ['FA', 'L1', 'L2', 'L3']:
    
            f_name= '/neurospin/grip/external_databases/dHCP_CR_JD_2018/release3/dhcp_dmri_shard_pipeline/sub-{}/ses-{}/dwi/DTI/dtifit_b1000/sub-{}_ses-{}_{}.nii.gz'.format(
                                    subj_id, ses_id,subj_id, ses_id, metric)

            vol = aims.read(f_name)
            im = vol.arraydata()[0]
    
            for vox in uniques:
                im[vox[0],vox[1],vox[2]] = 0
        
            #new_im = merged_vol.astype(np.int16)


            new_im = aims.Volume(im)
            new_im.header().update(vol.header())

            aims.write(new_im, f_name)

            if os.path.isfile(f_name):
                print('{} {} {} metric finished'.format(subj_id, ses_id, metric))
    else: 
        print('{} {} NO ABERRANT VALUES FOUND'.format(subj_id, ses_id))
    print('***')

### ===================== RUN over SUBJECTS =============== ##
if __name__ == "__main__":
	print('START at: {}'.format(datetime.datetime.now()))

	### read in the file containing the list of subject ID and session IDs
	if len(sys.argv) < 2:
		print("You must provide subject file!")
		sys.exit()
	else:
		subject_file = sys.argv[1]

	subjects = pd.read_csv(subject_file, names=['subject_id', 'session_id', 'template'])
	
	for i, row in subjects.iterrows():
		
		subject_id = row.subject_id
		session_id = row.session_id
		
		print('WORKING ON SUBJECT {}, ID {}'.format(subject_id, session_id))


		iDIR='/neurospin/grip/external_databases/dHCP_CR_JD_2018/release3/dhcp_dmri_shard_pipeline/sub-{}/ses-{}/dwi/DTI/dtifit_b1000'.format(subject_id, session_id)

		print(os.path.join(iDIR, 'sub-{}_ses-{}_{}.nii.gz'.format(subject_id, session_id, 'FA')))
		## check DTI exists
		if os.path.exists(os.path.join(iDIR, 'sub-{}_ses-{}_{}.nii.gz'.format(subject_id, session_id, 'FA'))):
			
			correct_DTI(subj_id=subject_id, ses_id=session_id)
		
		else:
        		print('{} {} NO DTI... Skipping'.format(subject_id, session_id))
        		print('***')
	print('END at: {}'.format(datetime.datetime.now()))





