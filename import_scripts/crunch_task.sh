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
BLOCK=$3

if test $SUB -eq "010"; then
	RUN=$((RUN+1))
fi
	
# CRUNCH Specific Paths
DICOMDIR_A=/mnt/hgfs/ORIGDATA/Day_A/
DICOMDIR_B=/mnt/hgfs/ORIGDATA/Day_B/

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


for dir in ${DICOMDIR_B}/*
do
	dir=${dir%*/}
	dir=${dir##*/}

	if [ "${dir/$SUB}" != "$dir" ] ; then
		DICOMDIR_B=$DICOMDIR_B/$dir
	fi
done


	
{
for F in "${DICOMDIR_A}/*TFE*.PAR"
do
	{
	 	python3.6 -W ignore /home/fsluser/Documents/rs_proc/import_scripts/par2nii.py -i $F -o ${OUTDIR}/t1.nii
	} &> /dev/null
done
} || {
	for F in "${DICOMDIR_A}/*TFE*.nii"
	do
		cp $F ${OUTDIR}/t1.nii
	done

}




{
	python3.6 -W ignore /home/fsluser/Documents/rs_proc/import_scripts/par2nii_multiblock.py -i ${DICOMDIR_B} -p 'exp_run' -t 'nii' -b ${BLOCK} -o ${OUTDIR}/func_data.nii

} || {
	python3.6 -W ignore /home/fsluser/Documents/rs_proc/import_scripts/par2nii_multiblock.py -i ${DICOMDIR_B} -p 'exp_run' -t 'PAR' -b ${BLOCK} -o ${OUTDIR}/func_data.nii
}

for F in "${DICOMDIR_B}/*fieldmap*.PAR"
do

	if test -f ${F}; then
		python3.6 -W ignore /home/fsluser/Documents/rs_proc/import_scripts/par2nii.py -i $F -o ${OUTDIR}/fmap.nii
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
