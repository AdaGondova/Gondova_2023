#!/bin/bash

### project FA maps to NeonatAtlas template space for the QC 
### requires fsl_init


while read -r line
do 
	
	subject_id=$(echo "${line}" | cut -d ',' -f1)
	session_id=$(echo "${line}" | cut -d ',' -f2)
	template_age=$(echo "${line}" | cut -d ',' -f3)

	### Reference
	### will depends on AGE
	iREF=/neurospin/grip/external_databases/dHCP_CR_JD_2018/Projects/eLife_replication/SourceData/Neonate_Atlas/Single_Subject/${template_age}wk_Single_Subject/${template_age}wk_DTI_SS_FA.nii 

	### native FA
	iDiffFolder=/neurospin/grip/external_databases/dHCP_CR_JD_2018/release3/dhcp_dmri_shard_pipeline/sub-${subject_id}/ses-${session_id}/dwi
	iFA=${iDiffFolder}/DTI/dtifit_b1000/sub-${subject_id}_ses-${session_id}_FA.nii.gz
	
	### native2template non-linear warp
	iWARP=/neurospin/grip/external_databases/dHCP_CR_JD_2018/Projects/eLife_replication/DerivedData/subjects/sub-${subject_id}/ses-${session_id}/sub-${subject_id}_ses-${session_id}_shard2NeonateA_nonlinear_transf.nii.gz

	### outdata in /Results 
	oFA=/neurospin/grip/external_databases/dHCP_CR_JD_2018/Projects/eLife_replication/Results/registration_QC/sub-${subject_id}_ses-${session_id}_FA_in_Natlas_template_${template_age}.nii.gz

	if [ ! -f "$oFA" ]
	then 
		echo Warping ${subject_id} ${session_id} ...
		applywarp --ref=${iREF} --in=${iFA} --warp=${iWARP} --out=${oFA}
	else
		echo ${subject_id} ${session_id} already registered to template.

	fi 

	if [ -f "$oFA" ]		
	then 
		echo DONE
	fi 
done <  "$1"
echo FINISHED AT: `date`

