#!/bin/bash

# requires fsl_init

echo Stared at: `date`

## subject info / line 
while read -r line
do
	subject=$(echo "${line}" | cut -d ',' -f1)
	session=$(echo "${line}" | cut -d ',' -f2)
	echo WORKING ON SUBJECT: "$subject", session "$session"



	### =================== REGISTER shard 2 anat ====== ###
	iREF=/neurospin/grip/external_databases/dHCP_CR_JD_2018/release3/dhcp_anat_pipeline/sub-${subject}/ses-${session}/anat/sub-${subject}_ses-${session}_desc-restore_T2w.nii.gz
	iWARP=/neurospin/grip/external_databases/dHCP_CR_JD_2018/release3/dhcp_dmri_shard_pipeline/sub-${subject}/ses-${session}/xfm/sub-${subject}_ses-${session}_from-dwi_to-T2w_mode-image.mat

	## will be created and crushed in the loop to save space
	mkdir -p tmp
	

	### REGISTER metrics 
	for metric in FA L1 #MD RD
	do 
		iMETRIC=/neurospin/grip/external_databases/dHCP_CR_JD_2018/release3/dhcp_dmri_shard_pipeline/sub-${subject}/ses-${session}/dwi/DTI/dtifit_b1000/sub-${subject}_ses-${session}_${metric}.nii.gz
		oMETRIC=tmp/registered_${metric}.nii.gz

		#registration 
		echo Registering ${metric}
		flirt -in ${iMETRIC} -ref ${iREF} -out ${oMETRIC} -init ${iWARP} -applyxfm -interp trilinear

		if [ -f "$oMETRIC" ]
		then 
			echo ${metric} REGISTERED
		fi	
	done 	
	
	### REGISTER cortical parcellation 
	#iSEGM=/neurospin/grip/external_databases/dHCP_CR_JD_2018/Projects/eLife_replication/DerivedData/subjects/sub-${subject}/ses-${session}/sub-${subject}_ses-${session}_NAt_cortex_shard_space.nii.gz
	#oSEGM=tmp/registered_cortex.nii.gz

	#if [ -f "$iSEGM" ]
	#then 	
	#	echo skipping
		#flirt -in ${iSEGM} -ref ${iREF} -out ${oSEGM} -init ${iWARP} -applyxfm -interp nearestneighbour
	#fi 

	#if [ -f "$oSEGM" ]
	#then 
			#echo ${oSEGM} REGISTERED

			## treat segmentations 
			### create hemi masks, then voronoi
			
			### left 
			#echo creating left hemi mask
			#AimsThreshold -i /neurospin/grip/external_databases/dHCP_CR_JD_2018/release3/dhcp_anat_pipeline/sub-${subject}/ses-${session}/anat/sub-${subject}_ses-${session}_desc-ribbon_dseg.nii.gz -m be -t 2 -u 3 -o tmp/mask_L.nii.gz -b
			#AimsMask -i tmp/registered_cortex.nii.gz -o tmp/masked_cortex_L.nii.gz -m tmp/mask_L.nii.gz
			#echo computing left voronoi
			#AimsFileConvert -i tmp/masked_cortex_L.nii.gz -o tmp/masked_cortex_L_S16.nii.gz -t S16
			#AimsVoronoi -i tmp/masked_cortex_L_S16.nii.gz -o tmp/registered_cortex_voronoi_L.nii.gz

			### right
			#echo creating right hemi mask
			#AimsThreshold -i /neurospin/grip/external_databases/dHCP_CR_JD_2018/release3/dhcp_anat_pipeline/sub-${subject}/ses-${session}/anat/sub-${subject}_ses-${session}_desc-ribbon_dseg.nii.gz -m be -t 41 -u 42 -o tmp/mask_R.nii.gz -b
			#AimsMask -i tmp/registered_cortex.nii.gz -o tmp/masked_cortex_R.nii.gz -m tmp/mask_R.nii.gz
			#echo computing right voronoi
			#AimsFileConvert -i tmp/masked_cortex_R.nii.gz -o tmp/masked_cortex_R_S16.nii.gz -t S16
			#AimsVoronoi -i tmp/masked_cortex_R_S16.nii.gz -o tmp/registered_cortex_voronoi_R.nii.gz

	#fi


	### =============================== EXTRACTION ========== ###

	OUTDIR=/neurospin/grip/external_databases/dHCP_CR_JD_2018/Projects/eLife_replication/DerivedData/subjects/sub-${subject}/ses-${session}
	
	for hemi in left right 
	do 
		if [ "${hemi}" == left ]
		then 
			in_hemi=L
		else
			in_hemi=R
		fi 

		#for metric in cortex_voronoi #FA RD MD L1
		#for metric in FA RD MD L1
		for metric in FA
		do 
			echo $hemi $metric ...
			iMESH=/neurospin/grip/external_databases/dHCP_CR_JD_2018/release3/dhcp_anat_pipeline/sub-${subject}/ses-${session}/anat/sub-${subject}_ses-${session}_T2w_${in_hemi}white_bv_transformed.gii
			#oTEXTURE=${OUTDIR}/sub-${subject}_ses-${session}_${hemi}_texture_${metric}_majority.gii
			oTEXTURE=${OUTDIR}/sub-${subject}_ses-${session}_${hemi}_texture_${metric}.gii

			if [ "${metric}" == cortex_voronoi ]
			then 
				metric=cortex_voronoi_${in_hemi}
			fi
			echo ${metric}
			
			##/i2bm/brainvisa/brainvisa-master/bin/bv AimsVol2Tex -i ${oMETRIC} -m ${iMESH} -o ${oTEXTURE} -a ${iANOTHER} -height 1.5 -radius 1.5 -v 3
			# OK one
			#/i2bm/brainvisa/brainvisa-master/bin/bv AimsVol2Tex -a tmp/registered_${metric}.nii.gz -m ${iMESH} -o ${oTEXTURE} -i tmp/registered_${metric}.nii.gz -height 1.5  -radius 0.75 -v 3				
			/i2bm/brainvisa/brainvisa-master/bin/bv AimsVol2Tex -a tmp/registered_${metric}.nii.gz -m ${iMESH} -o ${oTEXTURE} -i tmp/registered_L1.nii.gz -height 1.5 -radius 1.5 -v 3
			
			if [ -f "$oTEXTURE" ]
			then 
				echo $oTEXTURE created
			fi 
		done 
	done 
	echo "$subject" "$session" DONE
	rm -r tmp
done < "$1"
echo Finished at: `date`





