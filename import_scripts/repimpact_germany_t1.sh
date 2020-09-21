#!/bin/bash

# Copies necessary files to processing [OUTDIR] directory
# !!!IMPORTANT!!!
#
# This is what OUTDIR is expected to contain after importing
#
# t1.nii  	 -> Highres anatomical file
# func_data.nii  -> 4D EPI file
#
# [OPTIONAL]
#
# rev_phase.nii  -> 4D EPI acquired with reverse phase [for fieldmap correction]
SUB=$1
OUTDIR=$2
DATADIR=/home/luna.kuleuven.be/u0101486/workspace/data/RepImpact/Raw/GermanyNIFTIfiles/$SUB

if [ ! -d "$DATADIR" ]; then
	echo "Cannot find ${SUB}. Data import failed."
fi





	for F in "${DATADIR}/*T1W*.nii.gz"
	do
		cp $F ${OUTDIR}/t1_${SUB}.nii.gz
		gunzip -f ${OUTDIR}/t1_${SUB}.nii.gz
	done
