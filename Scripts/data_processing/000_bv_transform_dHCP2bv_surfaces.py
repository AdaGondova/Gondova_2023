import os, datetime, subprocess, sys, io
import pandas as pd 
import numpy as np 
from soma import aims

### this needs debugging!

### script to transform dhcp surfaces to work with AIMS 
### not pretty but works
### requires bv bash, 
### expects subject file

def get_trm_file(subject_id, session_id):
	
	iDIR= '/neurospin/grip/external_databases/dHCP_CR_JD_2018/release3/dhcp_anat_pipeline/sub-{}/ses-{}/anat'.format(subject_id, session_id)

	## files to convert
	t2w = os.path.join(iDIR, 'sub-{}_ses-{}_desc-restore_T2w.nii.gz'.format(subject_id, session_id))	
	t2w_new = os.path.join(iDIR, 'sub-{}_ses-{}_T2w_copy.nii.gz'.format(subject_id, session_id))
	conv_out = os.path.join(iDIR, 'sub-{}_ses-{}_T2w.ima'.format(subject_id, session_id))


	## first needs converting to older image type to get transform info (!!! to keep the original)
	command = subprocess.Popen(['cp', t2w, t2w_new], shell=False)
	command.wait()
	command = subprocess.Popen(['AimsFileConvert',
                             '-i',
                              t2w_new,
                             '-o',
                              conv_out], shell=False)
	command.wait()

	if os.path.isfile(conv_out+'.minf'):
		file = io.open(conv_out+'.minf', 'r', encoding = "ISO-8859-1")
	else:
		print('{} {} info file not created, something is wrong!'.format(subject_id, session_id))
		return
	
	
	for line in file.readlines():
		if 'transformations' in line:
			line = line.strip()
			line = line[line.find('[')+3 : line.find(']')]
			trm =  np.asarray([np.float64(x) for x in line.split(',')]).reshape(4,4)[:-1]
			trm = np.insert(trm[:,:-1], 0, trm[:,-1], 0) 
    	
			to_save = os.path.join(iDIR, 'sub-{}_ses-{}.trm'.format(subject_id, session_id))
			np.savetxt(to_save, trm, delimiter=" ", fmt='%1.16f')
            
			if os.path.isfile(to_save):
				print('TRM SAVED')
			else: 
				print('Something went wrong.')
 
	file.close()
	return to_save

def transform_surfaces(subject_id, session_id, path_to_trm):

	iDIR = '/neurospin/grip/external_databases/dHCP_CR_JD_2018/release3/dhcp_anat_pipeline/sub-{}/ses-{}/anat'.format(subject_id, session_id)

	Rwhite = os.path.join(iDIR, 'sub-{}_ses-{}_hemi-right_wm.surf.gii'.format(subject_id, session_id))
	Lwhite = os.path.join(iDIR, 'sub-{}_ses-{}_hemi-left_wm.surf.gii'.format(subject_id, session_id))

	to_save = os.path.join(iDIR, 'sub-{}_ses-{}.trm'.format(subject_id, session_id))
	
	surfaces = [Rwhite, Lwhite]
	out_names = [ os.path.join(iDIR, 'sub-{}_ses-{}_T2w_Rwhite_bv_transformed.gii'.format(subject_id, session_id)),
			os.path.join(iDIR, 'sub-{}_ses-{}_T2w_Lwhite_bv_transformed.gii'.format(subject_id, session_id)) ]  

	for file, output in zip(surfaces, out_names):
		print('Transforming {} ...'.format(file))

		if not os.path.exists(output):
			command = subprocess.Popen(['AimsMeshTransform',
                             '-i',
                              file,
                             '-o',
                             output,
                             '-t',
                            to_save ], shell=False)
			command.wait()
		else:
			print('{} already exists'.format(output))
	
	if os.path.isfile(out_names[0]):
		print('{} {} DONE'.format(subject_id, session_id))
		return True
	else:
		print('Transformation did not work! INVESTIGATE')
		return False


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
		# transform surfaces
		path_to_trm = get_trm_file(subject_id, session_id)
		print(path_to_trm)
		## this is strange!
		if path_to_trm is not None:
			transform = transform_surfaces(subject_id=subject_id, session_id=session_id, path_to_trm=path_to_trm)

	
	print('END at: {}'.format(datetime.datetime.now()))	




