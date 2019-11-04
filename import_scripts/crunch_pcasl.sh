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
DICOMDIR_A=/mnt/hgfs/ORIGDATA/Day_A/

if [ ! -d "${OUTDIR}" ]; then
	mkdir -m 777 ${OUTDIR}
else
	{
	rm -f -r ${OUTDIR}/*
	} &> /dev/null
fi



for dir in ${DICOMDIR_A}/*
do
	dir=${dir%*/}
	dir=${dir##*/}

	if [ "${dir/$SUB}" != "$dir" ] ; then
		DICOMDIR_A=$DICOMDIR_A/$dir
	fi
done


for F in "${DICOMDIR_A}/*pCASL*.PAR"
do
	echo $F
	#dcm2niix_afni -s y -v 2 -p n -o ${OUTDIR} $F
	#3dPAR2AFNI.pl -n -o ${OUTDIR} $F
	python par2nii.py -i $F -o ${OUTDIR}/pCASL.nii

done

	
{
for F in "${DICOMDIR_A}/*TFE*.PAR"
do
	{
		dcm2niix_afni -s y -v 2 -p n -o ${OUTDIR} $F
	} &> /dev/null
done
} || {
	for F in "${DICOMDIR_A}/*TFE*.nii"
	do
		cp $F ${OUTDIR}
	done

}



for F in "${OUTDIR}/*TFE*.nii"
do
	mv $F ${OUTDIR}/t1.nii
done


mkdir ${OUTDIR}/perf
mkdir ${OUTDIR}/perf/rest
mkdir ${OUTDIR}/perf/anat
mv ${OUTDIR}/pCASL.nii ${OUTDIR}/perf/rest/pCASL.nii
mv ${OUTDIR}/t1.nii ${OUTDIR}/perf/anat/t1_structural.nii





matlab -nodesktop -nosplash -softwareopengl <<<"addpath(genpath('/home/fsluser/Documents/spm12')); batch_run('${OUTDIR}');"


# Delete intermediate files
rm ${OUTDIR}/perf/rest/pCASL.nii
rm ${OUTDIR}/perf/rest/sASLflt_rpCASL.nii
rm ${OUTDIR}/perf/rest/ASLflt_rpCASL.nii
rm ${OUTDIR}/perf/rest/rpCASL.nii


