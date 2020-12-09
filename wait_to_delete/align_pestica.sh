#!/bin/bash

AFNI_BYTEORDER=LSB_FIRST
export AFNI_BYTEORDER

PDIR='/home/luna.kuleuven.be/u0101486/workspace/fmri_proc/ext/pestica_afni_v5.2/'
SDIR='/home/luna.kuleuven.be/u0101486/workspace/fmri_proc/ext/pestica_afni_v5.2/slomoco/'

OUTDIR=$1
FNAME=$2
TR=$3
DO_SLOMOCO=$4

rm -f ${OUTDIR}/nat_mask.nii
3dAutomask -prefix ${OUTDIR}/nat_mask.nii ${OUTDIR}/${FNAME}.nii

source ${PDIR}/setup_pestica.sh

cd $OUTDIR

rm -f -r ${OUTDIR}/PESTICA5
rm -f -r ${OUTDIR}/SLOMOCO5


cp -f ${OUTDIR}/slice_acq_s.txt ${OUTDIR}/tshiftfile_sec.1D
cp -f ${OUTDIR}/slice_acq_ms.txt ${OUTDIR}/tshiftfile.1D

cp -f ${OUTDIR}/slice_acq_ms.txt ${OUTDIR}/PESTICA5/tshiftfile.1D
cp -f ${OUTDIR}/slice_acq_s.txt ${OUTDIR}/PESTICA5/tshiftfile_sec.1D

${PDIR}/run_pestica.sh -d func_data.nii -s "1 2 3 4 5" -b -m 2
${SDIR}/slicemoco_newalgorithm.sh -d func_data.nii -r -m 2

rm -f ${OUTDIR}/r${FNAME}.nii
3dcopy ${OUTDIR}/SLOMOCO5/func_data.slicemocoxy_afni.slomoco.pestica+orig ${OUTDIR}/r${FNAME}.nii
copy_header ${OUTDIR}/${FNAME}.nii ${OUTDIR}/r${FNAME}.nii
