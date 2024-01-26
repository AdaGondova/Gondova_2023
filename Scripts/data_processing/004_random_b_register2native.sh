#!/bin/bash

### Registers random parcellations from dHCP template to native space

# to use wb_command
. /volatile/dhcp-structural-pipeline/parameters/path.sh 


echo Stared at: `date`

## Loop over subjects! (assumes subj and session IDs are in comma sep csv!) 
while read -r line
do
	subject=$(echo "${line}" | cut -d ',' -f1)
	session=$(echo "${line}" | cut -d ',' -f2)

	INPUTDIR=/neurospin/grip/external_databases/dHCP_CR_JD_2018/release3/dhcp_anat_pipeline/sub-${subject}/ses-${session}
	## SAVE IT TO eLIFE Derived Data 
	OUTDIR=/neurospin/grip/external_databases/dHCP_CR_JD_2018/Projects/eLife_replication/DerivedData/subjects/sub-${subject}/ses-${session}

	echo Working on: "$subject" "$session"


	for hemi in left right ; do
			echo Working on "$hemi"

			iNatSurf=${INPUTDIR}/anat/sub-${subject}_ses-${session}_hemi-${hemi}_wm.surf.gii
			iNatTransSphere=${INPUTDIR}/xfm/sub-${subject}_ses-${session}_hemi-${hemi}_from-native_to-dhcpSym40_dens-32k_mode-sphere.surf.gii
			
			echo INPUT: ${iNatSurf}

			iTemplateSurf=/neurospin/grip/external_databases/dHCP_CR_JD_2018/Projects/andrea/SourceData/dhcpSym_template/week-40_hemi-${hemi}_space-dhcpSym_dens-32k_wm.surf.gii
			iTemplateSphere=/neurospin/grip/external_databases/dHCP_CR_JD_2018/Projects/andrea/SourceData/dhcpSym_template/week-40_hemi-${hemi}_space-dhcpSym_dens-32k_sphere.surf.gii

			iNatMidthickness=${INPUTDIR}/anat/sub-${subject}_ses-${session}_hemi-${hemi}_midthickness.surf.gii
			iTemplateMidthickness=/neurospin/grip/external_databases/dHCP_CR_JD_2018/Projects/andrea/SourceData/dhcpSym_template/week-40_hemi-${hemi}_space-dhcpSym_dens-32k_midthickness.surf.gii

			#####change here 	
			for parc in 64 128 256 512 ; do 
				iTemplateParc=/neurospin/grip/external_databases/dHCP_CR_JD_2018/Projects/eLife_replication/DerivedData/random_parcellation/${hemi}_template_random_parc_pathKmeans_${parc}_clusters.tex.gii
				oNativeParc=${OUTDIR}/sub-${subject}_ses-${session}_hemi-${hemi}_pathKmeans_${parc}.label.gii

				#check files exist 
				check=true
				for iFile in $iNatSurf $iNatTransSphere $iTemplateSurf $iTemplateSphere $iNatMidthickness $iTemplateMidthickness $iTemplateParc
				do 
					if [ ! -e $iFile ]; then 
						echo ${iFile} not found
						check=false
					fi
				done

				if [ "$check"  = false ]; then 
					echo Input files not found. Skipping "$subject" "$session" "$hemi" parce: ${parc}
				else
					wb_command -label-resample "$iTemplateParc" "$iTemplateSphere" "$iNatTransSphere" ADAP_BARY_AREA "$oNativeParc" -area-surfs "$iTemplateMidthickness" "$iNatMidthickness"
				fi 

				if [ -e $oNativeParc ]; then 
					echo "$oNativeParc" written. 
				fi 
			done
	done	
done < "$1"
echo Finished at: `date`
