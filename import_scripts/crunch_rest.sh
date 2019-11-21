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
# TODO CORRIGIR!
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


for F in "${DICOMDIR_A}/*GRASE*.PAR"
do
	echo $F
	#dcm2niix_afni -s y -v 2 -p n -o ${OUTDIR} $F
	#3dPAR2AFNI.pl -n -o ${OUTDIR} $F
	python par2nii.py -i $F -o ${OUTDIR}/GRASE.nii

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


{
for F in ${DICOMDIR_A}/*RSN*.nii
do
	FNAME="$(basename -- $F)"
	RNUM=$(echo $FNAME| rev |cut -d'.' -f 2 | cut -d'_' -f 2 | rev )
	{
	cp $F ${OUTDIR} 
	} &> /dev/null
done
} || {
	for F in ${DICOMDIR_A}/*RSN*.PAR
	do
	{	
	dcm2niix_afni -s y -v 2 -p n -f '%s_%f_%p_%s' -o ${OUTDIR} $F
	} &> /dev/null
	done
}


for F in ${OUTDIR}/*RSN*.nii
do			
	mv $F ${OUTDIR}/func_data.nii
done


