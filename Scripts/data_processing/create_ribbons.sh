#!/bin/bash

#fsl_init
#. <>/dhcp-structural-pipeline-master/parameters/path.sh
#FSLDIR=/volatile/home/ag265252/fsl
#. ${FSLDIR}/etc/fslconf/fsl.sh
#PATH=${FSLDIR}/bin:${PATH}
#export FSLDIR PATH


LeftGreyRibbonValue="3"
LeftGreyRibbonValueIn="2"
RightGreyRibbonValue="42"
RightGreyRibbonValueIn="41"

CSF_label="1"
CGM_label="2"

### subject 

iDIR=/neurospin/grip/external_databases/dHCP_CR_JD_2018/release3/dhcp_anat_pipeline/sub-CC01145XX16/ses-98330/anat


T2=${iDIR}/sub-CC01145XX16_ses-98330_desc-restore_T2w.nii.gz
iLab=${iDIR}/sub-CC01145XX16_ses-98330_desc-drawem9_dseg.nii.gz
oRib=${iDIR}/sub-CC01145XX16_ses-98330_desc_ribbon_dseg.nii.gz



# create ribbon 
for h in left right ; do
	iPial=${iDIR}/sub-CC01145XX16_ses-98330_hemi-${h}_pial.surf.gii
	echo $iPial
	oDist=sub-CC01145XX16_ses-98330_dist_${h}.nii.gz

	#wb_command -create-signed-distance-volume $iPial $T2 $oDist
	echo calculated distance
	#fslmaths $oDist -uthr 0 -abs -bin $oDist	
	#fslmaths $iLab -mul $oDist sub-CC01145XX16_ses-98330_tissue_labels_${h}.nii.gz
	echo ${h} tissue labels created...
done	
## LH

L_tissue_labels=sub-CC01145XX16_ses-98330_tissue_labels_left.nii.gz

mirtk padding $L_tissue_labels $L_tissue_labels sub-CC01145XX16_ses-98330_in_left.nii.gz 2 $CSF_label $CGM_label 0
fslmaths sub-CC01145XX16_ses-98330_in_left.nii.gz -bin -mul $LeftGreyRibbonValueIn sub-CC01145XX16_ses-98330_in_left.nii.gz
fslmaths $L_tissue_labels -thr $CGM_label -uthr $CGM_label -bin -mul $LeftGreyRibbonValue sub-CC01145XX16_ses-98330_out_left.nii.gz

echo L tissue label created

## RH 

R_tissue_labels=sub-CC01145XX16_ses-98330_tissue_labels_right.nii.gz

mirtk padding $R_tissue_labels $R_tissue_labels sub-CC01145XX16_ses-98330_in_right.nii.gz 2 $CSF_label $CGM_label 0
fslmaths sub-CC01145XX16_ses-98330_in_right.nii.gz -bin -mul $RightGreyRibbonValueIn sub-CC01145XX16_ses-98330_in_right.nii.gz
fslmaths $R_tissue_labels -thr $CGM_label -uthr $CGM_label -bin -mul $RightGreyRibbonValue sub-CC01145XX16_ses-98330_out_right.nii.gz

echo R tissue label created

fslmaths sub-CC01145XX16_ses-98330_in_left.nii.gz -add sub-CC01145XX16_ses-98330_in_right.nii.gz -add sub-CC01145XX16_ses-98330_out_left.nii.gz -add sub-CC01145XX16_ses-98330_out_right.nii.gz $oRib

if [ -f "${oRib}" ]
then 
	echo Ribbon created 
fi 

