import numpy as np
import pandas as pd
from soma import aims
import os, datetime, subprocess, sys


### requires bv 

if __name__ == "__main__":
	print('START at: {}'.format(datetime.datetime.now()))

	labels = pd.read_csv('/neurospin/grip/external_databases/dHCP_CR_JD_2018/Projects/eLife_replication/SourceData/atlas_labelling_scheme.csv')
	sub_regions = labels['Label'].values

	
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
		print(subject_id)
		print(session_id)	

		iDIR = '/neurospin/grip/external_databases/dHCP_CR_JD_2018/Projects/eLife_replication/DerivedData/subjects/sub-{}/ses-{}'.format(subject_id, session_id)
		iParc = os.path.join(iDIR , 'sub-{}_ses-{}_NAt_shard_space.nii.gz'.format(subject_id, session_id))
		oParc = os.path.join(iDIR , 'sub-{}_ses-{}_NAt_cortex_shard_space.nii.gz'.format(subject_id, session_id))
		#print(iParc)

		if os.path.isfile(iParc):
			print('Extracting cortical parcellation...')
			parc = aims.read(iParc)
			parc_a = np.array(parc.arraydata(), dtype=int)

			parc_a[~np.isin(parc_a, sub_regions)] = 0
			parc_a = parc_a.astype(np.int16)

			oParcV = aims.Volume(parc_a)
			oParcV.header().update(parc.header())
			aims.write(oParcV, oParc)

		if os.path.isfile(oParc):
			print('Cortical parcellation created.')
	
	

	print('END at: {}'.format(datetime.datetime.now()))







