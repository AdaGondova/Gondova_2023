#!/bin/bash

#### run DTI on the shard dMRI
### inistialise FSL outside the script!
### accepts file with subject id, session id
# similar to eddy dHCP pipeline: https://git.fmrib.ox.ac.uk/matteob/dHCP_neo_dMRI_pipeline_release/-/blob/master/dHCP_neo_dMRI_runPostProc.sh 


while read -r line
do 
	
	subject_id=$(echo "${line}" | cut -d ',' -f1)
	session_id=$(echo "${line}" | cut -d ',' -f2)

	INPUTDIR=/neurospin/grip/external_databases/dHCP_CR_JD_2018/release3/dhcp_dmri_shard_pipeline/sub-${subject_id}/ses-${session_id}/dwi
	echo ${INPUTDIR}
	diffFolder=${INPUTDIR}/DTI

	if [ ! -d "$INPUTDIR" ]	
	then 
		echo ${subject_id} session ${session_id} NO DMRI!
	
	else

	   if [ -f "${diffFolder}/dtifit_b0/sub-${subject_id}_ses-${session_id}_FA.nii.gz" ]
	   then 

		echo ${subject_id} session ${session_id} already has DTI. CONTINUTE

	   else

		if [ ! -f "${INPUTDIR}/sub-${subject_id}_ses-${session_id}_desc-shard_motion.txt" ]
		then 
			echo ${subject_id} session ${session_id} does not have GOOD DTI
			echo ${subject_id},${session_id} >> dMRI_wrong_inputs_for_DTI.csv 


		else
	
			mkdir -p ${diffFolder}

			# create links to the input data
			ln -s ${INPUTDIR}/sub-${subject_id}_ses-${session_id}_desc-preproc_dwi.nii.gz ${diffFolder}/data.nii.gz
			ln -s ${INPUTDIR}/sub-${subject_id}_ses-${session_id}_desc-brain_mask.nii.gz ${diffFolder}/nodif_brain_mask.nii.gz
			ln -s ${INPUTDIR}/sub-${subject_id}_ses-${session_id}_desc-preproc_dwi.bval ${diffFolder}/bvals

			if [ -f "${INPUTDIR}/sub-${subject_id}_ses-${session_id}_desc-preproc_dwi_original.bvec" ] 
			then
				ln -s ${INPUTDIR}/sub-${subject_id}_ses-${session_id}_desc-preproc_dwi_original.bvec ${diffFolder}/bvecs
			else 
				ln -s ${INPUTDIR}/sub-${subject_id}_ses-${session_id}_desc-preproc_dwi.bvec ${diffFolder}/bvecs
			fi	

			echo Input creation finished
			echo Running DTI... 

			start_date=`date`

			
			#============================================================================
			# Run only on 0 and 1000b
			#
			#============================================================================
			uniqueBvals=(`cat ./shells`)
			#uniqueBvals=(`cat ${diffFolder}/bvals`)
			echo ${uniqueBvals}

			for b in "${uniqueBvals[@]}"
			do

			if [ $b == 0 ] || [ $b == 1000 ]
			then	
   	 			select_dwi_vols ${diffFolder}/data ${diffFolder}/bvals ${diffFolder}/mean_b${b} ${b} -m
   			
				echo "Fitting DT to b=${b} shell..."
		        	mkdir -p ${diffFolder}/dtifit_b${b}
	
		        	select_dwi_vols ${diffFolder}/data ${diffFolder}/bvals ${diffFolder}/dtifit_b${b}/b${b} 0 -b ${b} -obv ${diffFolder}/bvecs 
        			dtifit -k ${diffFolder}/dtifit_b${b}/b${b} -o ${diffFolder}/dtifit_b${b}/sub-${subject_id}_ses-${session_id} -m ${diffFolder}/nodif_brain_mask -r ${diffFolder}/dtifit_b${b}/b${b}.bvec -b ${diffFolder}/dtifit_b${b}/b${b}.bval --sse --save_tensor
        			fslmaths ${diffFolder}/mean_b${b} -div ${diffFolder}/mean_b0 -mul ${diffFolder}/nodif_brain_mask ${diffFolder}/att_b${b}
    				fslmaths ${diffFolder}/mean_b${b} -mul ${diffFolder}/nodif_brain_mask ${diffFolder}/mean_b${b}
			fi
			done
		fi
	  fi
	#============================================================================
	# Clean up files 
	#============================================================================
	echo cleaning up the folder...
	rm ${diffFolder}/att_b0.nii.gz
	rm ${diffFolder}/att_b1000.nii.gz
	rm ${diffFolder}/mean_b0.nii.gz
	rm ${diffFolder}/mean_b1000.nii.gz

	rm -r ${diffFolder}/dtifit_b0
	rm ${diffFolder}/dtifit_b1000/sub-${subject_id}_ses-${session_id}_V*
	rm ${diffFolder}/dtifit_b1000/sub-${subject_id}_ses-${session_id}_S*

	echo ${subject_id} finished...
	fi		
done <  "$1"
echo FINISHED AT: `date`
