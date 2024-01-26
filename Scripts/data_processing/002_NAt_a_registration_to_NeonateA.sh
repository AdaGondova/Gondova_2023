#!/bin/bash

### inistialise FSL outside the script!
### example here: https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FNIRT/UserGuide

##Registering FA-image to atlas
## FSL User Guide example:
#flirt -ref ${FSLDIR}/data/standard/FMRIB58_FA_1mm_brain -in my_FA -omat my_affine_transf.mat
#fnirt --in=my_FA --aff=my_affine_transf.mat --cout=my_nonlinear_transf --config=FA_2_FMRIB58_1mm
#applywarp --ref=${FSLDIR}/data/standard/FMRIB58_FA_1mm_brain --in=my_FA --warp=my_nonlinear_transf --out=my_warped_FA

#more /i2bm/local/fsl/etc/flirtsch/FA_2_FMRIB58_1mm.cnf

while read -r line
do 
	
	subject_id=$(echo "${line}" | cut -d ',' -f1)
	session_id=$(echo "${line}" | cut -d ',' -f2)
	template_age=$(echo "${line}" | cut -d ',' -f3)

	### will depends on AGE
	iREF=/neurospin/grip/external_databases/dHCP_CR_JD_2018/Projects/eLife_replication/SourceData/Neonate_Atlas/Single_Subject/${template_age}wk_Single_Subject/${template_age}wk_DTI_SS_FA.nii 
	iSEG=/neurospin/grip/external_databases/dHCP_CR_JD_2018/Projects/eLife_replication/SourceData/Neonate_Atlas/Single_Subject/${template_age}wk_Single_Subject/${template_age}wk_DTI_SS_atlas_labels.nii

	iDiffFolder=/neurospin/grip/external_databases/dHCP_CR_JD_2018/release3/dhcp_dmri_shard_pipeline/sub-${subject_id}/ses-${session_id}/dwi
	iFA=${iDiffFolder}/DTI/dtifit_b1000/sub-${subject_id}_ses-${session_id}_FA.nii.gz

	if [ ! -f "$iFA" ]
	then 
		echo ${subject_id} ${session_id} does not have FA data. SKIPPING...

	else

		#oDiffFolder=/neurospin/grip/external_databases/dHCP_CR_JD_2018/release3/dhcp_dmri_shard_pipeline/sub-${subject_id}/ses-${session_id}/dwi/NeonateA_reg

		oDiffFolder=/neurospin/grip/external_databases/dHCP_CR_JD_2018/Projects/eLife_replication/DerivedData/subjects/sub-${subject_id}/ses-${session_id}


		if [ ! -d "$oDiffFolder" ]
		then 
			mkdir -p ${oDiffFolder}
		fi 
		
		## Save to derived data instead???!!! 
		oAFF=${oDiffFolder}/sub-${subject_id}_ses-${session_id}_shard2NeonateA_affine_transf.mat
		oNL=${oDiffFolder}/sub-${subject_id}_ses-${session_id}_shard2NeonateA_nonlinear_transf
		oFinalNL=${oDiffFolder}/sub-${subject_id}_ses-${session_id}_NeonateA2shard_nonlinear_transf
		oFA=${oDiffFolder}/sub-${subject_id}_ses-${session_id}_NeonateA_FA.nii.gz
		oFA_flirt=${oDiffFolder}/sub-${subject_id}_ses-${session_id}_NeonateA_FA_flirt.nii.gz
		oParc=${oDiffFolder}/sub-${subject_id}_ses-${session_id}_NAt_shard_space.nii.gz

		if [ ! -f "${oParc}" ]
		then
			echo Running ${subject_id} ${session_id} registration...
			echo Running FLIRT...
			flirt -ref ${iREF} -in ${iFA} -omat ${oAFF} #-out ${oFA_flirt}
			echo Running FNIRT...
			fnirt --ref=${iREF} --in=${iFA} --aff=${oAFF} --cout=${oNL} --config=FA_2_FMRIB58_1mm_new ### original config overwrites reference images and makes everything broken!
		
			### test first without the inversion
			echo Inversing the warp 
			invwarp --ref=${iFA} --warp=${oNL} --out=${oFinalNL}
			echo Applying the inverse warp...
			applywarp --ref=${iFA} --in=${iSEG} --warp=${oFinalNL} --out=${oParc} --interp=nn

			if [ -f "$oParc" ]
			then 
				echo ${subject_id} ${session_id} DONE...
			fi	
		else 
			echo ${subject_id} ${session_id} already segmented. SKIPPING... 		
		fi

	fi	

	
	### register to the shard segmentation to anatomy
	#oDiffFolder=/neurospin/grip/external_databases/dHCP_CR_JD_2018/Projects/eLife_replication/DerivedData/segmentations/sub-${subject_id}/ses-${session_id}
	#iParc=${oDiffFolder}/sub-${subject_id}_ses-${session_id}_NAt_shard_space.nii.gz
	#oParc=${oDiffFolder}/sub-${subject_id}_ses-${session_id}_NAt_anat_space.nii.gz

	#iREF=/neurospin/grip/external_databases/dHCP_CR_JD_2018/release3/dhcp_anat_pipeline/sub-${subject}/ses-${session}/anat/sub-${subject}_ses-${session}_desc-restore_T2w.nii.gz
	#iWARP=/neurospin/grip/external_databases/dHCP_CR_JD_2018/release3/dhcp_dmri_shard_pipeline/sub-${subject}/ses-${session}/xfm/sub-${subject}_ses-${session}_from-dwi_to-T2w_mode-image.mat

	#if [ ! -f "$oParc" && -f "$iREF" ]
	#then 
	#	Registering to native...
	#	flirt -in ${iParc} -ref ${iREF} -out ${oParc} -init ${iWARP} -applyxfm -interp nn	
	#fi 
	

done <  "$1"
echo FINISHED AT: `date`



