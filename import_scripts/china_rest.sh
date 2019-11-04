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
	
	
# CRUNCH Specific Paths
DICOMDIR_T1=/mnt/hgfs/ssd_tmp/china/raw/t1/
DICOMDIR_RSN=/mnt/hgfs/ssd_tmp/china/raw/rest/

if [ ! -d "${OUTDIR}" ]; then
	mkdir -m 777 -p ${OUTDIR}
else
	{
	rm -f -r ${OUTDIR}/*
	} &> /dev/null
fi



for dir in ${DICOMDIR_T1}/*
do
	dir=${dir%*/}
	dir=${dir##*/}

	if [ "${dir/$SUB}" != "$dir" ] ; then
		DICOMDIR_T1=$DICOMDIR_T1/$dir
	fi
done



for dir in ${DICOMDIR_RSN}/*
do
	dir=${dir%*/}
	dir=${dir##*/}

	if [ "${dir/$SUB}" != "$dir" ] ; then
		DICOMDIR_RSN=$DICOMDIR_RSN/$dir
	fi
done
	

dcm2niix_afni -s n -v 2 -p n -o ${OUTDIR} ${DICOMDIR_T1}
dcm2niix_afni -s n -v 2 -p n -o ${OUTDIR} ${DICOMDIR_RSN}



for F in "${OUTDIR}/*t1_mprage*.nii"
do
	mv $F ${OUTDIR}/t1.nii
done

for F in "${OUTDIR}/*bold*.nii"
do
	mv $F ${OUTDIR}/func_data.nii
done

