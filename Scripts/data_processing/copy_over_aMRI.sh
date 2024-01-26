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

	if [ ! -f "$OUTDIR/anat/sub-${subject_id}_ses-${session_id}_desc-biasfield_T2w.nii.gz" ]
	
	then 	
		echo Copying $subject_id $session_id ... 

		mkdir -p /neurospin/grip/external_databases/dHCP_CR_JD_2018/release3/dhcp_anat_pipeline/sub-${subject_id}/ses-${session_id}/anat
		mkdir -p /neurospin/grip/external_databases/dHCP_CR_JD_2018/release3/dhcp_anat_pipeline/sub-${subject_id}/ses-${session_id}/xfm
		
		## copy t2 weighted
		cp ${INDIR}/anat/sub-${subject_id}_ses-${session_id}_desc-restore_T2w* ${OUTDIR}/anat/

		if [ ! -f "${INDIR}/anat/sub-${subject_id}_ses-${session_id}_hemi-left_wm.surf.gii" ]
		then
			echo $subject_id $session_id DOES NOT HAVE MESHES 
			echo $subject_id,$session_id >> missing_surface_data.csv

		else
			## copy meshes
			cp ${INDIR}/anat/sub-${subject_id}_ses-${session_id}_hemi-left_wm* ${OUTDIR}/anat/
			cp ${INDIR}/anat/sub-${subject_id}_ses-${session_id}_hemi-right_wm* ${OUTDIR}/anat/
			### needs midthickness!
			cp ${INDIR}/anat/sub-${subject_id}_ses-${session_id}_hemi-left_midthickness* ${OUTDIR}/anat/
			cp ${INDIR}/anat/sub-${subject_id}_ses-${session_id}_hemi-right_midthickness* ${OUTDIR}/anat/
		fi 
	
		## copy xfm data 
		cp ${INDIR}/xfm/* ${OUTDIR}/xfm/

		## copy info file
		cp /media/ag265252/JessicaDubois/rel3_dhcp_anat_pipeline/sub-${subject_id}/sub-${subject_id}_sessions.tsv /neurospin/grip/external_databases/dHCP_CR_JD_2018/release3/dhcp_anat_pipeline/sub-${subject_id}/
		
		if [ -f "$OUTDIR/anat/sub-${subject_id}_ses-${session_id}_desc-restore_T2w.nii.gz" ]
		then
			echo DONE
		fi
	else 
		echo $subject_id $session_id folder EXISTS!
	fi	
done <  "$1"

echo FINISHED AT: `date`


