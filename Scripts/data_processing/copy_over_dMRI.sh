#!/bin/bash

## copy over the shard data for the 379 subjects if not already on grip 

echo SUBJECTS TO COPY OVER:

while read -r line
do 
	
	subject_id=$(echo "${line}" | cut -d ',' -f1)
	session_id=$(echo "${line}" | cut -d ',' -f2)


	OUTDIR=/neurospin/grip/external_databases/dHCP_CR_JD_2018/release3/dhcp_dmri_shard_pipeline/sub-${subject_id}/ses-${session_id}
	INDIR=/media/ag265252/JessicaDubois/rel3_dhcp_dmri_shard_pipeline/sub-${subject_id}/ses-${session_id}/

	if [ ! -d "$OUTDIR" ]
	then 	
		echo Copying $subject_id $session_id ...

		mkdir -p /neurospin/grip/external_databases/dHCP_CR_JD_2018/release3/dhcp_dmri_shard_pipeline/sub-${subject_id}/	
		cp -R $INDIR $OUTDIR
		cp /media/ag265252/JessicaDubois/rel3_dhcp_dmri_shard_pipeline/sub-${subject_id}/sub-${subject_id}_sessions.tsv /neurospin/grip/external_databases/dHCP_CR_JD_2018/release3/dhcp_dmri_shard_pipeline/sub-${subject_id}/
		
		if [ -d "$OUTDIR/dwi" ]
		then
			echo DONE
		fi
	else 
		echo $subject_id $session_id folder EXISTS!
	fi	
done <  "$1"

echo FINISHED AT: `date`


