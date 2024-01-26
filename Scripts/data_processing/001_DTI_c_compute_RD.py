### Compute RD maps from the L2&L3 (after cleaning)

import numpy as np
import pandas as pd
from soma import aims
import os, datetime, subprocess, sys

### requires bv 

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
		print(subject_id)
		print(session_id)

		iDIR = iDIR = '/neurospin/grip/external_databases/dHCP_CR_JD_2018/release3/dhcp_dmri_shard_pipeline/sub-{}/ses-{}/dwi/DTI/dtifit_b1000'.format(subject_id, session_id)
		
		iL2 = os.path.join(iDIR, 'sub-{}_ses-{}_L2.nii.gz'.format(subject_id, session_id))
		iL3 = os.path.join(iDIR, 'sub-{}_ses-{}_L3.nii.gz'.format(subject_id, session_id))
		oRD = os.path.join(iDIR, 'sub-{}_ses-{}_RD.nii.gz'.format(subject_id, session_id))
		print(iL2)
		if os.path.isfile(iL2) and os.path.isfile(iL3):
			print('Computing RD')

			L2 = aims.read(iL2)
			L2_a = np.array(L2.arraydata())

			L3 = aims.read(iL3)
			L3_a = np.array(L3.arraydata())

			RD_a = (L2_a + L3_a)/2
			oRD_V = aims.Volume(RD_a)
			oRD_V.header().update(L2.header())
			aims.write(oRD_V, oRD)

		if os.path.isfile(oRD):
			print('RD computed.')

	print('END at: {}'.format(datetime.datetime.now()))
			
