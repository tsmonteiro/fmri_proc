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
DATADIR=/home/luna.kuleuven.be/u0101486/workspace/data/RepImpact/Raw/BelgiumNIFTIfiles/NIFTIfiles/$SUB

if [ ! -d "$DATADIR" ]; then
	echo "Cannot find ${SUB}. Data import failed."
else
	if [ ! -d "${OUTDIR}" ]; then
		mkdir -p -m 777 ${OUTDIR}
	else
		{
		rm ${OUTDIR}/*
		} &> /dev/null
	fi



	for F in "${DATADIR}/*rsn_2*.nii.gz"
	do
		cp $F ${OUTDIR}/func_data.nii.gz
		gunzip -f ${OUTDIR}/func_data.nii.gz
	done

	{
	for F in "${DATADIR}/*rsn_rev*.nii.gz"
	do
		cp $F ${OUTDIR}/rev_phase.nii.gz
		gunzip -f ${OUTDIR}/rev_phase.nii.gz
	done
	} | {
	for F in "${DATADIR}/*rsn_REV*.nii.gz"
	do
		cp $F ${OUTDIR}/rev_phase.nii.gz
		gunzip -f ${OUTDIR}/rev_phase.nii.gz
	done

	}

  FLIST=`ls ${DATADIR}/*T1W*.nii.gz`

	for F in ${FLIST}
	do
    echo "Copying ${F} to ${OUTDIR}"
		cp $F ${OUTDIR}/t1.nii.gz
		gunzip -f ${OUTDIR}/t1.nii.gz
	done


fi
