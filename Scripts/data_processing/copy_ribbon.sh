#!/bin/bash

## copy over the shard data for the 379 subjects if not already on grip 

## to register diffusion to anatomy: required files:
## /neurospin/grip/external_databases/dHCP_CR_JD_2018/release3/dhcp_anat_pipeline/sub-${subject}/ses-${session}/anat/sub-${subject}_ses-${session}_desc-restore_T2w.nii.gz
## WM meshed- transformed??
## trm are available with diffusion data

echo SUBJECTS TO COPY OVER:

while read -r line
do 
	
	subject_id=$(echo "${line}" | cut -d ',' -f1)
	session_id=$(echo "${line}" | cut -d ',' -f2)


	OUTDIR=/neurospin/grip/external_databases/dHCP_CR_JD_2018/release3/dhcp_anat_pipeline/sub-${subject_id}/ses-${session_id}
	INDIR=/media/ag265252/JessicaDubois/rel3_dhcp_anat_pipeline/sub-${subject_id}/ses-${session_id}



	if [ ! -f "${OUTDIR}/anat/sub-${subject_id}_ses-${session_id}_desc-ribbon_dseg.nii.gz" ]
	then 

		echo Need to copy ${subject_id},${session_id}

		if [ -f "${INDIR}/anat/sub-${subject_id}_ses-${session_id}_desc-ribbon_dseg.nii.gz" ]
		then 
			echo copying ${subject_id},${session_id}
			cp ${INDIR}/anat/sub-${subject_id}_ses-${session_id}_desc-ribbon* ${OUTDIR}/anat/
		else
			echo ${subject_id},${session_id} > no_ribbon.csv
		fi

	fi 

done <  "$1"

echo FINISHED AT: `date`

	

