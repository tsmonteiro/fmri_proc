#!/bin/bash

AFNI_BYTEORDER=LSB_FIRST
export AFNI_BYTEORDER 

PDIR='/home/fsluser/Documents/rs_proc/ext/pestica4/'
SDIR='/home/fsluser/Documents/rs_proc/ext/pestica4/slomoco/'

OUTDIR=$1
FNAME=$2
TR=$3
DO_SLOMOCO=$4


3dcopy ${OUTDIR}/nat_mask.nii ${OUTDIR}/afni_data.brain+orig
3dcopy ${OUTDIR}/${FNAME}.nii ${OUTDIR}/afni_data+orig

source ${PDIR}/setup_pestica.sh

3drefit  -TR ${TR} \
	-Tslices `cat ${OUTDIR}/slice_acq.txt` \
       -denote ${OUTDIR}/afni_data+orig

mkdir ${OUTDIR}/pestica4
cp ${OUTDIR}/afni_data.brain+orig* ${OUTDIR}/pestica4/

echo ${OUTDIR}
cd $OUTDIR

${PDIR}/run_pestica.sh -d afni_data -s "1 2 3 4" -b

if [ "$DO_SLOMOCO" -eq "1" ]; then
	${SDIR}/slicemoco_newalgorithm.sh -d afni_data -r
	
	3dcopy ${OUTDIR}/pestica4/coupling_ret_card_pestica+orig ${OUTDIR}/coupling_ret_card_pestica.nii
	3dcopy ${OUTDIR}/pestica4/coupling_ret_resp_pestica+orig ${OUTDIR}/coupling_ret_resp_pestica.nii
	#TODO The p file is only necessary for comparisons... rmeove later
	3dcopy ${OUTDIR}/pestica4/afni_data.retroicor_pestica+orig ${OUTDIR}/p${FNAME}.nii
	3dcopy ${OUTDIR}/slomoco4/afni_data.slicemocoxy_afni.slomoco_pestica+orig ${OUTDIR}/rp${FNAME}.nii
else
	3dcopy ${OUTDIR}/pestica4/coupling_ret_card_pestica+orig ${OUTDIR}/coupling_ret_card_pestica.nii
	3dcopy ${OUTDIR}/pestica4/coupling_ret_resp_pestica+orig ${OUTDIR}/coupling_ret_resp_pestica.nii
	3dcopy ${OUTDIR}/pestica4/afni_data.retroicor_pestica+orig ${OUTDIR}/p${FNAME}.nii
fi


