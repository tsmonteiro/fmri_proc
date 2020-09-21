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
DICOMDIR_A=/media/u0101486/MONTEIRO/Day_A/
P2NII_PATH=/home/luna.kuleuven.be/u0101486/workspace/fmri_proc/import_scripts/


if [ ! -d "${OUTDIR}" ]; then
	mkdir -m 777 ${OUTDIR}
else
	{
	#rm -f -r ${OUTDIR}/*
  rm -f -r ${OUTDIR}/*
  rm -f -r ${OUTDIR}/FIX
  rm -f -r ${OUTDIR}/can_melodic.ic
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


#for F in "${DICOMDIR_A}/*pCASL*.PAR"
#do
#	echo $F
#	#dcm2niix_afni -s y -v 2 -p n -o ${OUTDIR} $F
#	#3dPAR2AFNI.pl -n -o ${OUTDIR} $F
#	python par2nii.py -i $F -o ${OUTDIR}/pCASL.nii
#
#done
#
#
#for F in "${DICOMDIR_A}/*GRASE*.PAR"
#do
#	echo $F
#	#dcm2niix_afni -s y -v 2 -p n -o ${OUTDIR} $F
#	#3dPAR2AFNI.pl -n -o ${OUTDIR} $F
#	python par2nii.py -i $F -o ${OUTDIR}/GRASE.nii
#
#done


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


for F in "${DICOMDIR_A}/*fieldmap*.PAR"
do

	if test -f ${F}; then
		python -W ignore ${P2NII_PATH}/par2nii.py -i $F -o ${OUTDIR}/fmap.nii
		HAS_FMAP=1

		{
			3dTcat -prefix ${OUTDIR}/fmap_mag.nii ${OUTDIR}/fmap.nii[0]
			3dTcat -prefix ${OUTDIR}/fmap_phase.nii ${OUTDIR}/fmap.nii[1]

			# FMAP has to be in rad/s, but
			# philips scanner saves phase image in Hz (1/deltaTE) --> -250Hz <-> 250Hz
			fslmaths ${OUTDIR}/fmap_phase -div 250 -mul 3.141593  -mul 1000 -div 3 ${OUTDIR}/fmap_phase_rads

			3dAutomask -prefix ${OUTDIR}/fmap_mask.nii -apply_prefix ${OUTDIR}/fmap_mag_brain.nii ${OUTDIR}/fmap_mag.nii

			fslmaths ${OUTDIR}/fmap_phase_rads -mas ${OUTDIR}/fmap_mask.nii ${OUTDIR}/fmap_phase_rads_m


			# A bit of regularization
			fugue --loadfmap=${OUTDIR}/fmap_phase_rads_m -s 1 --despike --savefmap=${OUTDIR}/fmap_phase_rads_m
			mv ${OUTDIR}/fmap_phase_rads_m.nii.gz ${OUTDIR}/fmap_phase_rads.nii.gz
		} &> /dev/null

	else
		HAS_FMAP=0
		echo "WARNING No FIELDMAP PAR/REC to be imported"
		touch ${OUTDIR}/no_fmap.txt
	fi
done
