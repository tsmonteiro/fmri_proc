#!/bin/bash
# Clear terminal screen
#printf "\033c"

SUB=$1
CONFIG=$2 #/home/fsluser/Documents/rs_proc/conf_files/rep_impact_belgium.conf
#BLOCK=$3
DOWN=$3

#TODO Add this to the configuration file
WDIR='/home/luna.kuleuven.be/u0101486/workspace/fmri_proc/'
QCDIR='/home/luna.kuleuven.be/u0101486/workspace/fmri_proc/QualityControl/'
DICERPATH='/home/luna.kuleuven.be/u0101486/workspace/fmri_proc/ext/DiCER/'
C3DPATH='/home/luna.kuleuven.be/u0101486/Software/itk/ITK/c3d/bin/'
CMTKPATH='/home/luna.kuleuven.be/u0101486/Software/CMTK/cmtk-3.3.1p1/build/bin/'

# PARSE Configuration file
ABIN=$(awk -F\=  '/^ABIN/{print $2}' "${CONFIG}" | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
FIXBIN=$(awk -F\=  '/^FIXBIN/{print $2}' "${CONFIG}" | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
EXTRACT_SURF=$(awk -F\=  '/^EXTRACT_SURF/{print $2}' "${CONFIG}" | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
DO_FUNC_SURF=$(awk -F\=  '/^DO_FUNC_SURF/{print $2}' "${CONFIG}" | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
MNI_REF=$(awk -F \= '/^\<MNI_REF\>/{print $2}' "${CONFIG}" | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
MNI_REF_2mm=$(awk -F\=  '/^\<MNI_REF_2mm\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
TMPDIR=$(awk -F\=  '/^\<TMPDIR\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
OUTDIR=$(awk -F\=  '/^\<OUTDIR\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
FINALDIR=$(awk -F\=  '/^\<FINALDIR\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
DO_COPY=$(awk -F\=  '/^\<DO_COPY\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
MOCO=$(awk -F\=  '/^\<MOCO\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
ICA_TYPE=$(awk -F\=  '/^\<ICA_TYPE\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
DO_REG=$(awk -F\=  '/^\<DO_REG\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
DO_REG2=$(awk -F\=  '/^\<DO_REG2\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
DO_ICA=$(awk -F\=  '/^\<DO_ICA\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
DO_NORM=$(awk -F\=  '/^\<DO_NORM\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
DO_QA=$(awk -F\=  '/^\<DO_QA\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
DO_QA_NATIVE=$(awk -F\=  '/^\<DO_QA_NATIVE\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
DO_CLEAN=$(awk -F\=  '/^\<DO_CLEAN\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
APPLY_FIX=$(awk -F\=  '/^\<APPLY_FIX\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
APPLY_DICER=$(awk -F\=  '/^\<APPLY_DICER\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
DO_FUNC_BIAS=$(awk -F\=  '/^\<DO_FUNC_BIAS\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
DO_ATLAS_NAT=$(awk -F\=  '/^\<DO_ATLAS_NAT\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
TR=$(awk -F\=  '/^\<TR\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
ACQ_TYPE=$(awk -F\=  '/^\<ACQ_TYPE\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
F2A_FUNC=$(awk -F\=  '/^\<F2A_FUNC\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
DO_SLC=$(awk -F\=  '/^\<DO_SLC\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
DO_FMAP=$(awk -F\=  '/^\<DO_FMAP\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
EES=$(awk -F\=  '/^\<EES\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
DO_PEST=$(awk -F\=  '/^\<DO_PEST\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
DO_ANAT_FMAP=$(awk -F\=  '/^\<DO_ANAT_FMAP\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
CORRECT_ANAT_FMAP_NATIVE=$(awk -F\=  '/^\<CORRECT_ANAT_FMAP_NATIVE\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
DO_DEOBL=$(awk -F\=  '/^\<DO_DEOBL\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
DO_DPK=$(awk -F\=  '/^\<DO_DPK\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
DO_ING=$(awk -F\=  '/^\<DO_ING\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
GLOB_VAL=$(awk -F\=  '/^\<GLOB_VAL\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')

FIX_CLASSIFIER=$(awk -F\=  '/^\<FIX_CLASSIFIER\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
FIX_CL_LABEL=$(awk -F\=  '/^\<FIX_CL_LABEL\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
FIX_THR=$(awk -F\=  '/^\<FIX_THR\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')


BREXT_TYPE=$(awk -F\=  '/^\<BREXT_TYPE\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')

DO_FWHM=$(awk -F\=  '/^\<DO_FWHM\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
DO_NLM=$(awk -F\=  '/^\<DO_NLM\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
USE_TISSUE_PRIOR=$(awk -F\=  '/^\<USE_TISSUE_PRIOR\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
GM_PRIOR=$(awk -F\=  '/^\<GM_PRIOR\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
WM_PRIOR=$(awk -F\=  '/^\<WM_PRIOR\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
CSF_PRIOR=$(awk -F\=  '/^\<CSF_PRIOR\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
REG_MODEL=$(awk -F\=  '/^\<REG_MODEL\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
EXTRACT_NUIS=$(awk -F\=  '/^\<EXTRACT_NUIS\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
FD_THR=$(awk -F\=  '/^\<FD_THR\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
IMPORT_SCRIPT=$(awk -F\=  '/^\<IMPORT_SCRIPT\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')


FINALDIR=`echo ${FINALDIR//@SUB@/${SUB}} | tr -d '"'`
OUTDIR=`echo ${OUTDIR//@SUB@/${SUB}} | tr -d '"'`
TMPDIR=`echo ${TMPDIR} | tr -d '"'`
ABIN=`echo ${ABIN} | tr -d '"'`
ANAT_PROC=`echo ${ANAT_PROC} | tr -d '"'`




# In case data has multiple blocks (multiple resting state blocks or perhaps task)
#if test ! -z "${BLOCK}"; then
#      FINALDIR=${FINALDIR}_${BLOCK}
#      OUTDIR=${OUTDIR}_${BLOCK}
#fi

#${CMTKPATH}/asegment ${OUTDIR}/anat_ref.nii atlas_image.nii atlas_labels.nii output_labels.ni
#exit

# Run the command below to train the classifier
#${FIXBIN}/fix -t /home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/fix_classifier_rs/classifier -l  ${OUTDIR}/
# exit

print_debug()
{
      MSG=$1
      MSG_CMD=$2
      HOUR=`date +%T`
      DAY=`date +%F`
      # This makes easier to find the message in log files
      echo ''
      echo '//////////////////////////////////////////////////////////////////////////////////'
      echo '//////////////////////////////////////////////////////////////////////////////////'
      echo ''

      echo "[DEBUG] (${DAY} - ${HOUR} ) ${MSG}"
      echo "$MSG_CMD"

      echo ''
      echo '//////////////////////////////////////////////////////////////////////////////////'
      echo '//////////////////////////////////////////////////////////////////////////////////'
      echo ''
      echo ''
}


log_command_div()
{
      HOUR=`date +%T`
      DAY=`date +%F`
      MSG=$1

      echo ''
      echo ''
      echo ''
      echo '*****************************************************************************************'
      echo '*****************************************************************************************'
      echo ''
      echo "(${DAY} - ${HOUR} ) $1 DONE"
      echo ''
      echo '*****************************************************************************************'
      echo '*****************************************************************************************'
      echo '*****************************************************************************************'
      echo ''
      echo ''

}

estimate_func_bias()
{

      TN=$1
      OUTDIR=$2
      ABIN=$3

      if (( $TN < 10 )); then
            FI=00$TN
      elif (( $TN < 100 )); then
            FI=0$TN
      else
            FI=$TN
      fi

      {
      ${ABIN}/N4BiasFieldCorrection -i ${OUTDIR}/tmp/func_data.${FI}.nii -s 3  -o [ ${OUTDIR}/tmp/bcorr_${FI}.nii,${OUTDIR}/tmp/biasfield_${FI}.nii ]
      } &> /dev/null
}
export -f estimate_func_bias

nlm_smooth()
{

      TN=$1
      OUTDIR=$2
      ABIN=$3
      DO_NLM=$4



      if (( $TN < 10 )); then
            FI=00$TN
      elif (( $TN < 100 )); then
            FI=0$TN
      else
            FI=$TN
      fi

      {

          ${ABIN}/DenoiseImage -x ${OUTDIR}/mask_mni.nii -n Gaussian -i ${OUTDIR}/tmp/func_data.${FI}.nii  -o [ ${OUTDIR}/tmp/func_data_n.${FI}.nii ]
      } &> /dev/null

}
export -f nlm_smooth

NVOLS=`fslnvols ${OUTDIR}/proc_data_native_fix.nii`

GDIR=/home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/StructAna

#rm -f ${OUTDIR}/t1_head.nii
#rm -f ${OUTDIR}/t1_head_nlm.nii
#rm -r -f ${OUTDIR}/mri
#rm -r -f ${OUTDIR}/report
#rm -r -f ${OUTDIR}/label

#ANATOMY processing
ANATDIR=${OUTDIR}/anat
#mkdir ${ANATDIR}
rm -f ${ANATDIR}/*

cp ${OUTDIR}/t1.nii ${ANATDIR}/anat.nii
cp ${OUTDIR}/proc_data_native.nii ${ANATDIR}/func_data.nii

#FUNC
fslmaths ${ANATDIR}/func_data.nii -Tmean ${ANATDIR}/mean_func_data.nii
gunzip -f ${ANATDIR}/mean_func_data.nii.gz
${ABIN}/DenoiseImage -i ${ANATDIR}/mean_func_data.nii -r 4x4x4 -n Rician -o [ ${ANATDIR}/mean_func_data_n.nii ]
3dUnifize -EPI -prefix ${ANATDIR}/mean_func_data_nd.nii ${ANATDIR}/mean_func_data_n.nii

#ANAT
3dresample -prefix ${ANATDIR}/anat_orig.nii -orient lpi -input ${ANATDIR}/anat.nii
${ABIN}/N4BiasFieldCorrection -i ${ANATDIR}/anat_orig.nii -s 3  -o [ ${ANATDIR}/anat_wb.nii ]
python ${WDIR}/extract_brain.py -i "${ANATDIR}/anat_wb.nii" -o "${ANATDIR}/brain_mask.nii"
${ABIN}/DenoiseImage -i ${ANATDIR}/anat_wb.nii -r 4x4x4 -n Rician -x ${ANATDIR}/brain_mask.nii -o [ ${ANATDIR}/anat_wbn.nii ]
3dcalc -a ${ANATDIR}/anat_wbn.nii -b ${ANATDIR}/brain_mask.nii -expr 'a*b' -prefix ${ANATDIR}/anat_proc_brain.nii
3dUnifize -prefix ${ANATDIR}/anat_ref.nii ${ANATDIR}/anat_wbn.nii
mv ${ANATDIR}/anat_wbn.nii ${ANATDIR}/anat_proc.nii
rm -f ${ANATDIR}/anat_w*.nii
rm -f ${ANATDIR}/anat.nii
exit
#cp ${ANATDIR}/anat_proc.nii /home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/Template/anat/${SUB}.nii


exit

#mkdir ${OUTDIR}/subseg
#rm -f ${OUTDIR}/subseg/*

#rm -f ${OUTDIR}/subseg/proc_data_mni_tmp.nii
#3dTproject -prefix ${OUTDIR}/subseg/proc_data_mni_tmp.nii -tr ${TR} -polort 1 -stopband 0 0.009  \
#          -ort ${OUTDIR}/acompcor_fix.txt \
#          -input ${OUTDIR}/proc_data_mni.nii

#3dmaskave -nball 10 -16 8 8 ${OUTDIR}/subseg/proc_data_mni_tmp.nii |  awk '{print $1}' > ${OUTDIR}/subseg/r_thalamus.txt
#3dmaskave -nball -10 -16 8 8 ${OUTDIR}/subseg/proc_data_mni_tmp.nii |  awk '{print $1}' > ${OUTDIR}/subseg/l_thalamus.txt
#3dmaskave -q -nball -22 8 4 8 ${OUTDIR}/subseg/proc_data_mni_tmp.nii > ${OUTDIR}/subseg/l_putamen.txt
#3dmaskave -q -nball 26 8 4 8 ${OUTDIR}/subseg/proc_data_mni_tmp.nii > ${OUTDIR}/subseg/r_putamen.txt
#3dmaskave -q -nball 0 -66 38 8 ${OUTDIR}/subseg/proc_data_mni_tmp.nii > ${OUTDIR}/subseg/pcc.txt


#3dTcorr1D -prefix ${OUTDIR}/subseg/Corr_L_Thal.nii -Fisher -mask ${OUTDIR}/mask_mni.nii ${OUTDIR}/subseg/proc_data_mni_tmp.nii ${OUTDIR}/subseg/l_thalamus.txt
#3dTcorr1D -prefix ${OUTDIR}/subseg/Corr_R_Thal.nii -Fisher -mask ${OUTDIR}/mask_mni.nii ${OUTDIR}/subseg/proc_data_mni_tmp.nii ${OUTDIR}/subseg/r_thalamus.txt
#3dTcorr1D -prefix ${OUTDIR}/subseg/Corr_R_Puta.nii -Fisher -mask ${OUTDIR}/mask_mni.nii ${OUTDIR}/subseg/proc_data_mni_tmp.nii ${OUTDIR}/subseg/l_putamen.txt
#3dTcorr1D -prefix ${OUTDIR}/subseg/Corr_L_Puta.nii -Fisher -mask ${OUTDIR}/mask_mni.nii ${OUTDIR}/subseg/proc_data_mni_tmp.nii ${OUTDIR}/subseg/r_putamen.txt
#3dTcorr1D -prefix ${OUTDIR}/subseg/Corr_PCC.nii -Fisher -mask ${OUTDIR}/mask_mni.nii ${OUTDIR}/subseg/proc_data_mni_tmp.nii ${OUTDIR}/subseg/pcc.txt
#rm -f /home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/VBM/s${SUB}.nii

#SIGMA=2.55 #6mm
#SIGMA=3.45 #8mm
#fslmaths ${OUTDIR}/subseg/Corr_L_Thal.nii   -kernel gauss ${SIGMA} -fmean /home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/VBM/s${SUB}.nii
#gunzip -f /home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/VBM/s${SUB}.nii.gz

#rm -f ${OUTDIR}/subseg/proc_data_mni_tmp.nii
#exit


#rm -f /home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/VBM/s${SUB}.nii


rm -f ${OUTDIR}/proc_data_mni_tmp.nii
#tcompcor_fix
#fslmaths ${OUTDIR}/mri/mwp1t1_head.nii -thr 0.1 -bin ${OUTDIR}/mri/gm_mask_mni
3dTproject -prefix ${OUTDIR}/proc_data_mni_tmp.nii -tr ${TR} -polort 1  -passband 0.01 0.08 -norm \
          -ort ${OUTDIR}/tcompcor_fix.txt \
          -input ${OUTDIR}/proc_data_mni.nii
rm -f /home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/REHO_New/${SUB}.nii
rm -f ${OUTDIR}/reho_tmp.nii
#${OUTDIR}/reho_tmp.nii
#-mask ${OUTDIR}/mask_mni.nii
3dReHo -inset ${OUTDIR}/proc_data_mni_tmp.nii  -prefix  ${OUTDIR}/reho_tmp.nii #/home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/REHO_New/${SUB}.nii
rm -f ${OUTDIR}/proc_data_mni_tmp.nii

#wm_sig.txt

#SIGMA=2.55 #6mm
SIGMA=3.45 #8mm
fslmaths ${OUTDIR}/reho_tmp.nii   -kernel gauss ${SIGMA} -fmean /home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/REHO_New/${SUB}.nii
gunzip -f /home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/REHO_New/${SUB}.nii.gz
#3dmerge -prefix /home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/REHO/aall_${SUB}.nii -1blur_fwhm 6 ${OUTDIR}/mri/reho_cat_tmp.nii
exit

#fslmaths  ${OUTDIR}/proc_data_mni.nii -Tmean ${OUTDIR}/mri/mean_func_mni.nii
#
#3dAllineate -input ${OUTDIR}/mri/mean_func_mni.nii.gz -base ${OUTDIR}/mri/wmt1_head.nii -1Dmatrix_save ${OUTDIR}/mri/reho_mat \
#    -prefix ${OUTDIR}/mni2anat.nii -onepass -cmass -cost lpc+hel
#
#
3dAllineate -input ${OUTDIR}/mri/reho_cat_tmp.nii -base ${OUTDIR}/mri/wmt1_head.nii -1Dmatrix_apply ${OUTDIR}/mri/reho_mat.aff12.1D \
        -prefix ${OUTDIR}/mri/reho_cat_tmp2.nii -onepass -cost ls

3dcalc -a ${OUTDIR}/mri/reho_cat_tmp2.nii -b ${OUTDIR}/mri/gm_mask_mni.nii.gz -c  ${OUTDIR}/mri/wj_t1_head.nii -expr '(a*b)' -prefix ${OUTDIR}/mri/reho_cat_all_anat.nii
#
#rm -f ${OUTDIR}/mri/reho_cat_tmp.nii
#rm -f ${OUTDIR}/mri/reho_cat_tmp2.nii
#rm -f ${OUTDIR}/proc_data_mni_tmp.nii
#
#
#3dTproject -prefix ${OUTDIR}/proc_data_mni_tmp.nii -tr ${TR} -polort 2 -passband 0.05 0.12 -norm \
#          -ort ${OUTDIR}/acompcor_fix.txt \
#          -input ${OUTDIR}/proc_data_mni.nii

#3dReHo -inset ${OUTDIR}/proc_data_mni_tmp.nii -prefix ${OUTDIR}/mri/reho_cat_tmp.nii

#fslmaths  ${OUTDIR}/proc_data_mni.nii -Tmean ${OUTDIR}/mri/mean_func_mni.nii

#3dAllineate -input ${OUTDIR}/mri/mean_func_mni.nii.gz -base ${OUTDIR}/mri/wmt1_head.nii -1Dmatrix_save ${OUTDIR}/mri/reho_mat \
#    -prefix ${OUTDIR}/mni2anat.nii -onepass -cmass -cost lpc+hel
#
#
#3dAllineate -input ${OUTDIR}/mri/reho_cat_tmp.nii -base ${OUTDIR}/mri/wmt1_head.nii -1Dmatrix_apply ${OUTDIR}/mri/reho_mat.aff12.1D \
#        -prefix ${OUTDIR}/mri/reho_cat_tmp2.nii -onepass -cmass -cost lpc+hel
#
#3dcalc -a ${OUTDIR}/mri/reho_cat_tmp2.nii -b ${OUTDIR}/mri/gm_mask_mni.nii.gz -expr 'a*b' -prefix ${OUTDIR}/mri/reho_cat_hi_anat.nii

#rm -f ${OUTDIR}/mri/reho_cat_tmp.nii
#rm -f ${OUTDIR}/mri/reho_cat_tmp2.nii
#rm -f ${OUTDIR}/proc_data_mni_tmp.nii
#

#rm -f /home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/REHO/${SUB}_ahi.nii
#rm -f /home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/REHO/${SUB}_alo.nii
#rm -f /home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/REHO/ahi_${SUB}.nii
#rm -f /home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/REHO/alo_${SUB}.nii
rm -f /home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/REHO/aall_${SUB}.nii
3dmerge -prefix /home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/REHO/aall_${SUB}.nii -1blur_fwhm 6 ${OUTDIR}/mri/reho_cat_all_anat.nii
rm -f ${OUTDIR}/mri/reho_cat_all_anat.nii
#3dmerge -prefix /home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/REHO/alo_${SUB}.nii -1blur_fwhm 8 ${OUTDIR}/mri/reho_cat_lo_anat.nii

#cp -f  ${OUTDIR}/mri/reho_cat_hi_anat.nii /home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/REHO/${SUB}_ahi.nii
#cp -f  ${OUTDIR}/mri/reho_cat_lo_anat.nii /home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/REHO/${SUB}_alo.nii

exit
fslmaths ${OUTDIR}/mri/mwp1t1_head.nii -thr 0.1 -bin ${OUTDIR}/mri/gm_mask_mni
3dReHo -inset ${OUTDIR}/proc_data_mni.nii -prefix ${OUTDIR}/reho_cat.nii -mask ${OUTDIR}/mri/gm_mask_mni.nii.gz
exit


fslmaths ${OUTDIR}/mri/wj_t1_head.nii -kernel gauss 2.5 -fmean /home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/VBM/s${SUB}.nii
gunzip -f /home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/VBM/s${SUB}.nii.gz
MET=fast
fslmaths ${OUTDIR}/mri/p0t1_head.nii -thr 1.8 -uthr 2.2 -bin ${OUTDIR}/mri/gm_mask.nii
fslmaths ${OUTDIR}/first/SubcortSeg_all_${MET}_firstseg.nii.gz -mas ${OUTDIR}/mri/wm_mask.nii.gz ${OUTDIR}/first/tmp_subcort.nii

fslstats ${OUTDIR}/first/tmp_subcort.nii -l 9.5 -u 11.5 -V |  awk '{print $2}' > ${GDIR}/L_Thal_Vol_${SUB}_CAT_WM.txt
fslstats ${OUTDIR}/first/tmp_subcort.nii -l 48.5 -u 50.5 -V |  awk '{print $2}' > ${GDIR}/R_Thal_Vol_${SUB}_CAT_WM.txt
fslstats ${OUTDIR}/first/tmp_subcort.nii -l 11.5 -u 13.5 -V |  awk '{print $2}' > ${GDIR}/L_Puta_Vol_${SUB}_CAT_WM.txt
fslstats ${OUTDIR}/first/tmp_subcort.nii -l 50.5 -u 52.5 -V |  awk '{print $2}' > ${GDIR}/R_Puta_Vol_${SUB}_CAT_WM.txt


exit

#3dresample -prefix ${OUTDIR}/t1_head.nii -orient lpi -input ${OUTDIR}/t1.nii
UNLMDIR=/home/luna.kuleuven.be/u0101486/workspace/matlab/toolbox/UnbiasedNonLocalMeans/bin/Linux/
${UNLMDIR}/UnbiasedNonLocalMeans --sigma 10 --rs 4,4,4  ${OUTDIR}/t1_head.nii ${OUTDIR}/t1_head_nlm.nii
#matlab "-nodesktop -nosplash " <<<"cd '${WDIR}'; cat12_seg('${OUTDIR}', '${OUTDIR}/t1_head.nii'); exit;"
#exit
rm -f ${OUTDIR}/first/tmp_subcort.nii
rm -f ${OUTDIR}/first/tmp_subcort.nii.gz
echo ${SUB}
fslmaths ${OUTDIR}/mri/p0t1_head.nii -thr 1.66 -uthr 2.5 -bin ${OUTDIR}/mri/gm_mask.nii
fslmaths ${OUTDIR}/mri/p0t1_head.nii -thr 2.8 -uthr 3.01 -bin ${OUTDIR}/mri/wm_mask.nii
fslmaths ${OUTDIR}/mri/p0t1_head.nii -thr 0.5 -uthr 1.2 -bin ${OUTDIR}/mri/csf_mask.nii
MET=fast
rm -f ${OUTDIR}/t1_tmp.nii
3dcalc -a ${OUTDIR}/t1_head_nlm.nii -b ${OUTDIR}/AnatMask.nii -expr 'a*b' -prefix ${OUTDIR}/t1_tmp.nii
#-a ${OUTDIR}/anat2atlas_fsl.mat
run_first_all -m ${MET}  -a ${OUTDIR}/anat2atlas_fsl.mat -i ${OUTDIR}/t1_tmp.nii -o ${OUTDIR}/first/SubcortSeg



#3dcalc -a ${OUTDIR}/first/SubcortSeg_all_${MET}_firstseg.nii.gz -b ${OUTDIR}/mri/gm_mask.nii.gz -expr 'a*b' -prefix ${OUTDIR}/first/tmp_subcort.nii

#fslmaths ${OUTDIR}/first/SubcortSeg_all_${MET}_firstseg.nii.gz -mas ${OUTDIR}/gm_mask.nii ${OUTDIR}/first/tmp_subcort.nii
cp -f ${OUTDIR}/first/SubcortSeg_all_${MET}_firstseg.nii.gz ${OUTDIR}/first/tmp_subcort.nii.gz

#cp ${OUTDIR}/mri/mwp1t1_head.nii /home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/VBM/${SUB}.nii
rm -f /home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/VBM/s${SUB}.nii
#3dmerge -prefix /home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/VBM/s${SUB}.nii -1blur_fwhm 7 ${OUTDIR}/mri/mwp1t1_head.nii
fslmaths ${OUTDIR}/mri/mwp1t1_head.nii -kernel gauss 3 -fmean /home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/VBM/s${SUB}.nii
gunzip -f /home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/VBM/s${SUB}.nii.gz

fslstats ${OUTDIR}/mri/p0t1_head.nii -V |  awk '{print $2}' > /home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/VBM/${SUB}_TIV.txt

fslstats ${OUTDIR}/mri/p0t1_head.nii -V |  awk '{print $2}' > ${GDIR}/TIV_${SUB}_CAT.txt
fslstats ${OUTDIR}/gm_mask.nii -V |  awk '{print $2}' > ${GDIR}/GM_${SUB}_CAT.txt
fslstats ${OUTDIR}/wm_mask.nii -V |  awk '{print $2}' > ${GDIR}/WM_${SUB}_CAT.txt
fslstats ${OUTDIR}/csf_mask.nii -V |  awk '{print $2}' > ${GDIR}/CSF_${SUB}_CAT.txt

fslstats ${OUTDIR}/first/tmp_subcort.nii -l 9.5 -u 11.5 -V |  awk '{print $2}' > ${GDIR}/L_Thal_Vol_${SUB}_CAT.txt
fslstats ${OUTDIR}/first/tmp_subcort.nii -l 48.5 -u 50.5 -V |  awk '{print $2}' > ${GDIR}/R_Thal_Vol_${SUB}_CAT.txt
fslstats ${OUTDIR}/first/tmp_subcort.nii -l 11.5 -u 13.5 -V |  awk '{print $2}' > ${GDIR}/L_Puta_Vol_${SUB}_CAT.txt
fslstats ${OUTDIR}/first/tmp_subcort.nii -l 50.5 -u 52.5 -V |  awk '{print $2}' > ${GDIR}/R_Puta_Vol_${SUB}_CAT.txt
rm -f ${OUTDIR}/first/tmp_subcort.nii.gz

cp ${OUTDIR}/report/cat_t1_head.mat ${GDIR}/${SUB}_CAT.mat
cp ${OUTDIR}/label/catROI_t1_head.mat ${GDIR}/${SUB}_ROI_CAT.mat

rm -f ${GDIR}/T1_${SUB}.mat
cp  ${OUTDIR}/t1_tmp.nii ${GDIR}/T1_${SUB}.nii

cp ${OUTDIR}/first/SubcortSeg_all_${MET}_firstseg.nii.gz ${GDIR}/${SUB}_FirstSeg.nii.gz

cp ${OUTDIR}/first/SubcortSeg-L_Puta_first.bvars ${GDIR}/L_Puta_${SUB}.bvars
cp ${OUTDIR}/first/SubcortSeg-R_Puta_first.bvars ${GDIR}/R_Puta_${SUB}.bvars
cp ${OUTDIR}/first/SubcortSeg-L_Thal_first.bvars ${GDIR}/L_Thal_${SUB}.bvars
cp ${OUTDIR}/first/SubcortSeg-R_Thal_first.bvars ${GDIR}/R_Thal_${SUB}.bvars


cp ${OUTDIR}/first/SubcortSeg-L_Puta_first.vtk ${GDIR}/L_Puta_${SUB}.vtk
cp ${OUTDIR}/first/SubcortSeg-R_Puta_first.vtk ${GDIR}/R_Puta_${SUB}.vtk
cp ${OUTDIR}/first/SubcortSeg-L_Thal_first.vtk ${GDIR}/L_Thal_${SUB}.vtk
cp ${OUTDIR}/first/SubcortSeg-R_Thal_first.vtk ${GDIR}/R_Thal_${SUB}.vtk



rm -f ${OUTDIR}/t1_tmp.nii

rm -f ${OUTDIR}/first/tmp_subcort.nii
exit
ROI_LABELS=${OUTDIR}/label/catROI_t1_head_nlm.xml


VGM=0
VWM=0
while IFS= read -r line
do
  if [ "${VGM}" -eq "0" ] || [ "${VWM}" -eq "0" ]; then
        if [[ $line == *"<Vgm>"* ]] && [ "${VGM}" -eq "0" ]; then
          echo "$line" > ${GDIR}/All_Vol_GM_${SUB}_CAT.txt
          VGM=1
        fi

        if [[ $line == *"<Vgm>"* ]] && [ "${VGM}" -eq "1" ]; then
          echo "$line" > ${GDIR}/All_Vol_GMN_${SUB}_CAT.txt
          VGM=1
        fi

        if [[ $line == *"<Vwm>"* ]]; then
          echo "$line" > ${GDIR}/All_Vol_WM_${SUB}_CAT.txt
          VWM=1
        fi
  fi


done < "$ROI_LABELS"



#REF=${WDIR}/atlas/MNI152_T1_1mm_brain.nii

#mkdir ${OUTDIR}/first_smooth

#rm -f -r ${OUTDIR}/first_smooth/
rm -f -r ${OUTDIR}/first/

mkdir ${OUTDIR}/first
cd ${OUTDIR}/first




#${ABIN}/N4BiasFieldCorrection -i ${OUTDIR}/t1_f.nii -s 4  -o [ ${OUTDIR}/t1_fb.nii,${OUTDIR}/biasfield.nii ]
#${C3DPATH}/c3d_affine_tool -ref ${REF} -src ${OUTDIR}/anat_ref.nii \
#                    -itk ${OUTDIR}/anat2atlas0GenericAffine.mat -ras2fsl -o ${OUTDIR}/anat2atlas_fsl.mat
#3dUnifize -prefix ${OUTDIR}/tmp_anat.nii ${OUTDIR}/anat_first.nii
UNLMDIR=/home/luna.kuleuven.be/u0101486/workspace/matlab/toolbox/UnbiasedNonLocalMeans/bin/Linux/

#unring.a64 ${OUTDIR}/t1_brain.nii ${OUTDIR}/tmp_anat.nii
#gunzip -f ${OUTDIR}/tmp_anat.nii.gz

#${UNLMDIR}/UnbiasedNonLocalMeans --sigma 13  ${OUTDIR}/t1_brain.nii ${OUTDIR}/anat_first.nii
#rm -f ${OUTDIR}/tmp_anat.nii
#fast --nobias -g  ${OUTDIR}/anat_first.nii

rm -f ${OUTDIR}/anat_first_gm_tmp.nii
rm -f ${OUTDIR}/anat_first_gm.nii.gz
rm -f ${OUTDIR}/anat_first_gm.nii

#3dcalc -a ${OUTDIR}/anat_first.nii -b ${OUTDIR}/anat_first_seg_1.nii.gz -expr 'a*b' -prefix ${OUTDIR}/anat_first_gm_tmp.nii
#fslmaths ${OUTDIR}/anat_first.nii -kernel 2D -fmedian ${OUTDIR}/anat_first_gm.nii.gz

run_first_all -v -b -m fast -s L_Puta,R_Puta,R_Thal,L_Thal -a ${OUTDIR}/anat2atlas_fsl.mat -i ${OUTDIR}/anat_first.nii -o ${OUTDIR}/first/SubcortSeg
#rm -f ${OUTDIR}/anat_first_gm.nii.gz

#exit
#3dUnifize -GM -prefix ${OUTDIR}/tmp_anat.nii ${OUTDIR}/anat_first.nii
#mv ${OUTDIR}/tmp_anat.nii ${OUTDIR}/anat_first.nii

#run_first_all -b -m fast -s L_Puta,R_Puta,R_Thal,L_Thal -a ${OUTDIR}/anat2atlas_fsl.mat -i ${OUTDIR}/anat_first.nii -o ${OUTDIR}/first/SubcortSeg

3dcalc -a ${OUTDIR}/first/SubcortSeg_all_fast_firstseg.nii.gz -b ${OUTDIR}/gm_mask.nii -expr 'a*b' -prefix ${OUTDIR}/first/tmp_subcort.nii


fslstats ${OUTDIR}/first/tmp_subcort.nii -l 9.5 -u 11.5 -V |  awk '{print $2}' > ${GDIR}/L_Thal_Vol_${SUB}.txt
fslstats ${OUTDIR}/first/tmp_subcort.nii -l 48.5 -u 50.5 -V |  awk '{print $2}' > ${GDIR}/R_Thal_Vol_${SUB}.txt
fslstats ${OUTDIR}/first/tmp_subcort.nii -l 11.5 -u 13.5 -V |  awk '{print $2}' > ${GDIR}/L_Puta_Vol_${SUB}.txt
fslstats ${OUTDIR}/first/tmp_subcort.nii -l 50.5 -u 52.5 -V |  awk '{print $2}' > ${GDIR}/R_Puta_Vol_${SUB}.txt
rm -f ${OUTDIR}/first/tmp_subcort.nii
exit

rm -f ${OUTDIR}/tmp_anat.nii

#run_first_all -b -m fast -s L_Puta,R_Puta,R_Thal,L_Thal -a ${OUTDIR}/anat2atlas_fsl.mat -i ${OUTDIR}/anat_ref.nii -o ${OUTDIR}/first_smooth/SubcortSeg

# Calculate volume

#10 Left-Thalamus-Proper 40
#11 Left-Caudate 30
#12 Left-Putamen 40
#13 Left-Pallidum 40
#16 Brain-Stem /4th Ventricle 40
#17 Left-Hippocampus 30
#18 Left-Amygdala 50
#26 Left-Accumbens-area 50
#49 Right-Thalamus-Proper 40
#50 Right-Caudate 30
#51 Right-Putamen 40
#52 Right-Pallidum 40
#53 Right-Hippocampus 30
#54 Right-Amygdala 50
#58 Right-Accumbens-area 50
#rm -f ${OUTDIR}/tmp.nii
fslstats ${OUTDIR}/AnatMask.nii -V |  awk '{print $2}' > ${GDIR}/TIV_${SUB}.txt
fslstats ${OUTDIR}/csf_mask.nii -V |  awk '{print $2}' > ${GDIR}/CSF_${SUB}.txt
fslstats ${OUTDIR}/wm_mask.nii -V |  awk '{print $2}' > ${GDIR}/WM_${SUB}.txt
fslstats ${OUTDIR}/gm_mask.nii -V |  awk '{print $2}' > ${GDIR}/GM_${SUB}.txt

fslstats ${OUTDIR}/first/SubcortSeg_all_fast_firstseg.nii.gz -l 9.5 -u 11.5 -V |  awk '{print $2}' > ${GDIR}/L_Thal_Vol_${SUB}.txt
fslstats ${OUTDIR}/first/SubcortSeg_all_fast_firstseg.nii.gz -l 48.5 -u 50.5 -V |  awk '{print $2}' > ${GDIR}/R_Thal_Vol_${SUB}.txt
fslstats ${OUTDIR}/first/SubcortSeg_all_fast_firstseg.nii.gz -l 11.5 -u 13.5 -V |  awk '{print $2}' > ${GDIR}/L_Puta_Vol_${SUB}.txt
fslstats ${OUTDIR}/first/SubcortSeg_all_fast_firstseg.nii.gz -l 50.5 -u 52.5 -V |  awk '{print $2}' > ${GDIR}/R_Puta_Vol_${SUB}.txt

cp ${OUTDIR}/first/SubcortSeg_all_fast_firstseg.nii.gz ${GDIR}/${SUB}_FirstSeg.nii.gz

cp ${OUTDIR}/first/SubcortSeg-L_Puta_first.bvars ${GDIR}/L_Puta_${SUB}.bvars
cp ${OUTDIR}/first/SubcortSeg-R_Puta_first.bvars ${GDIR}/R_Puta_${SUB}.bvars
cp ${OUTDIR}/first/SubcortSeg-L_Thal_first.bvars ${GDIR}/L_Thal_${SUB}.bvars
cp ${OUTDIR}/first/SubcortSeg-R_Thal_first.bvars ${GDIR}/R_Thal_${SUB}.bvars


cp ${OUTDIR}/first/SubcortSeg-L_Puta_first.vtk ${GDIR}/L_Puta_${SUB}.vtk
cp ${OUTDIR}/first/SubcortSeg-R_Puta_first.vtk ${GDIR}/R_Puta_${SUB}.vtk
cp ${OUTDIR}/first/SubcortSeg-L_Thal_first.vtk ${GDIR}/L_Thal_${SUB}.vtk
cp ${OUTDIR}/first/SubcortSeg-R_Thal_first.vtk ${GDIR}/R_Thal_${SUB}.vtk
exit
#cp ${OUTDIR}/t1_frn.nii /home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/StructAna/T1_${SUB}.nii
cp ${OUTDIR}/first/SubcortSeg-L_Puta_first.nii.gz /home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/StructAna/L_Puta_${SUB}.nii.gz
cp ${OUTDIR}/first/SubcortSeg-R_Puta_first.nii.gz /home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/StructAna/R_Puta_${SUB}.nii.gz
cp ${OUTDIR}/first/SubcortSeg-L_Thal_first.nii.gz /home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/StructAna/L_Thal_${SUB}.nii.gz
cp ${OUTDIR}/first/SubcortSeg-R_Thal_first.nii.gz /home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/StructAna/R_Thal_${SUB}.nii.gz


cp ${OUTDIR}/first/SubcortSeg-L_Puta_first.bvars /home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/StructAna/L_Puta_${SUB}.bvars
cp ${OUTDIR}/first/SubcortSeg-R_Puta_first.bvars /home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/StructAna/R_Puta_${SUB}.bvars
cp ${OUTDIR}/first/SubcortSeg-L_Thal_first.bvars /home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/StructAna/L_Thal_${SUB}.bvars
cp ${OUTDIR}/first/SubcortSeg-R_Thal_first.bvars /home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/StructAna/R_Thal_${SUB}.bvars


cp ${OUTDIR}/first/SubcortSeg-L_Puta_first.vtk /home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/StructAna/L_Puta_${SUB}.vtk
cp ${OUTDIR}/first/SubcortSeg-R_Puta_first.vtk /home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/StructAna/R_Puta_${SUB}.vtk
cp ${OUTDIR}/first/SubcortSeg-L_Thal_first.vtk /home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/StructAna/L_Thal_${SUB}.vtk
cp ${OUTDIR}/first/SubcortSeg-R_Thal_first.vtk /home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/StructAna/R_Thal_${SUB}.vtk
exit

#rm ${OUTDIR}/proc_data_mni_nlm.nii
NVOLS=`fslnvols ${OUTDIR}/proc_data_mni.nii`

#3dAutomask -prefix ${OUTDIR}/mask_mni.nii ${OUTDIR}/proc_data_mni.nii

#mkdir ${OUTDIR}/tmp
#3dTsplit4D -prefix ${OUTDIR}/tmp/func_data.nii -keep_datum ${OUTDIR}/proc_data_mni.nii
#parallel -j4 --line-buffer nlm_smooth ::: $(seq 0 ${NVOLS}) ::: ${OUTDIR} ::: ${ABIN}
#rm -f ${OUTDIR}/proc_data_mni_nlm.nii
#3dTcat -prefix ${OUTDIR}/proc_data_mni_nlm.nii -tr ${TR} ${OUTDIR}/tmp/func_data_n*.nii

rm -f ${OUTDIR}/gm_mask_mni_hi.nii
rm -f  ${OUTDIR}/gm_mask_mni_hi_tmp.nii.gz
rm -f  ${OUTDIR}/gm_mask_mni_hi_tmp.nii

rm -f ${OUTDIR}/gm_mask_mni_lo.nii
rm -f  ${OUTDIR}/gm_mask_mni_lo_tmp.nii.gz
rm -f  ${OUTDIR}/gm_mask_mni_lo_tmp.nii

rm -f ${OUTDIR}/csf_mask_mni_lo.nii
rm -f  ${OUTDIR}/csf_mask_mni_lo_tmp.nii.gz
rm -f  ${OUTDIR}/csf_mask_mni_lo_tmp.nii

fslmaths  ${OUTDIR}/gm_mask_mni.nii -thr 0.98 -bin -kernel 2D -dilM  ${OUTDIR}/gm_mask_mni_lo.nii
gunzip -f  ${OUTDIR}/gm_mask_mni_lo.nii.gz

fslmaths  ${OUTDIR}/csf_mask_mni.nii -thr 0.9 -bin -kernel 2D -ero  ${OUTDIR}/csf_mask_mni_lo.nii
gunzip -f  ${OUTDIR}/csf_mask_mni_lo.nii.gz

rm -r -f ${OUTDIR}/tmp
rm -f ${OUTDIR}/proc_data_mni_tmp.nii
rm -f  ${OUTDIR}/proc_data_mni_tmp_2.nii


python run_acompcor.py -d ${OUTDIR} -i proc_data_mni.nii -n csf_mask_mni_lo.nii -b mask_mni.nii -t ${TR} \
                      -aout acompcor_fix -tout tcompcor_fix -var 0.66

# Remove headers from compcor results
sed '1d' ${OUTDIR}/acompcor_fix.txt > ${OUTDIR}/tmp_acompcor_fix.txt
mv ${OUTDIR}/tmp_acompcor_fix.txt ${OUTDIR}/acompcor_fix.txt

sed '1d' ${OUTDIR}/tcompcor_fix.txt > ${OUTDIR}/tmp_tcompcor_fix.txt
mv ${OUTDIR}/tmp_tcompcor_fix.txt ${OUTDIR}/tcompcor_fix.txt


3dTproject -prefix ${OUTDIR}/proc_data_mni_tmp.nii -tr ${TR} -polort 2 -passband 0.01 0.15 -norm \
          -ort ${OUTDIR}/acompcor_fix.txt \
          -input ${OUTDIR}/proc_data_mni.nii

#3dTproject -prefix ${OUTDIR}/proc_data_mni_tmp_nlm.nii -tr ${TR} -polort 2 -stopband 0 0.01 -norm \
#          -ort ${OUTDIR}/tcompcor_fix.txt \
#          -input ${OUTDIR}/proc_data_mni.nii



##-neigh_RAD 3
3dReHo -inset ${OUTDIR}/proc_data_mni_tmp.nii -prefix ${OUTDIR}/reho.nii -mask ${OUTDIR}/gm_mask_mni_lo.nii
#3dReHo -inset ${OUTDIR}/proc_data_mni_tmp.nii -prefix ${OUTDIR}/reho_nlm.nii -neigh_RAD 3 -mask ${OUTDIR}/gm_mask_mni_lo.nii
rm -f ${OUTDIR}/reho_hi*.nii
rm -f ${OUTDIR}/reho_lo.nii
rm -f ${OUTDIR}/reho_lo_nlm.nii

3dcalc -a ${OUTDIR}/reho.nii -b ${OUTDIR}/gm_mask_mni_lo.nii -expr 'a*b' -prefix ${OUTDIR}/reho_lo.nii
#3dcalc -a ${OUTDIR}/reho_nlm.nii -b ${OUTDIR}/gm_mask_mni_lo.nii -expr 'a*b' -prefix ${OUTDIR}/reho_lo_nlm.nii

rm -f  ${OUTDIR}/reho.nii

#${ABIN}/DenoiseImage -x ${OUTDIR}/gm_mask_mni_hi.nii -n Gaussian -i ${OUTDIR}/reho_hi.nii -v 0 -o [ ${OUTDIR}/reho_hi_nlm.nii,  ${OUTDIR}/reho_hi_noise.nii ]

rm -f ${OUTDIR}/reho_lo_s*.nii

3dmerge -prefix ${OUTDIR}/reho_lo_s6.nii -1blur_fwhm 6 ${OUTDIR}/reho_lo.nii
rm -f  ${OUTDIR}/reho_1.nii
rm -f  ${OUTDIR}/reho_2.nii
rm -f  ${OUTDIR}/reho_3.nii

rm -f  ${OUTDIR}/proc_data_mni_tmp*.nii

exit

## L Putamen
#3dmaskave -q -nball -22 8 4 6 ${OUTDIR}/proc_data_mni.nii > ${OUTDIR}/l_putamen.txt
## R Putamen
#3dmaskave -q -nball 26 8 4 6 ${OUTDIR}/proc_data_mni.nii > ${OUTDIR}/r_putamen.txt
## R Putamen
#3dmaskave -q -nball 10 -16 8 6 ${OUTDIR}/proc_data_mni.nii > ${OUTDIR}/r_thalamus.txt
#3dmaskave -q -nball -16 -26 8 6 ${OUTDIR}/proc_data_mni.nii > ${OUTDIR}/l_thalamus.txt
#3dmaskave -q -nball -2 -50 -16 6 ${OUTDIR}/proc_data_mni.nii > ${OUTDIR}/crb_I_IV.txt
#rm -f ${OUTDIR}/Corr_L_Putamen.nii
#rm -f ${OUTDIR}/Corr_R_Putamen.nii
#rm -f ${OUTDIR}/Corr_L_Thalamus.nii
#rm -f ${OUTDIR}/Corr_R_Thalamus.nii
#rm -f ${OUTDIR}/Corr_CRB_I_V.nii
#3dTcorr1D -prefix ${OUTDIR}/Corr_L_Putamen.nii -mask ${OUTDIR}/mni_mask.nii.gz ${OUTDIR}/proc_data_mni.nii ${OUTDIR}/l_putamen.txt
#3dTcorr1D -prefix ${OUTDIR}/Corr_R_Putamen.nii -mask ${OUTDIR}/mni_mask.nii.gz ${OUTDIR}/proc_data_mni.nii ${OUTDIR}/r_putamen.txt
#3dTcorr1D -prefix ${OUTDIR}/Corr_L_Thalamus.nii -mask ${OUTDIR}/mni_mask.nii.gz ${OUTDIR}/proc_data_mni.nii ${OUTDIR}/l_thalamus.txt
#3dTcorr1D -prefix ${OUTDIR}/Corr_R_Thalamus.nii -mask ${OUTDIR}/mni_mask.nii.gz ${OUTDIR}/proc_data_mni.nii ${OUTDIR}/r_thalamus.txt
#3dTcorr1D -prefix ${OUTDIR}/Corr_CRB_I_V.nii -mask ${OUTDIR}/mni_mask.nii.gz ${OUTDIR}/proc_data_mni.nii ${OUTDIR}/crb_I_IV.txt
#rm -f  ${OUTDIR}/reho_1.nii
#rm -f  ${OUTDIR}/reho_2.nii
#rm -f  ${OUTDIR}/reho_3.nii
#exit


# =================================================
#
#            PREPROCESSING START
#
# =================================================
#if [ ! -d "${OUTDIR}" ]; then
#	mkdir -m 777 -p ${OUTDIR}
# mkdir -m 777 -p ${OUTDIR}/logs
# mkdir -m 777 -p ${OUTDIR}/images
#else
#	{
#	rm -f -r ${OUTDIR}/*
#	} &> /dev/null
#fi



if [ "$DO_COPY" -eq "1" ]; then
      printf " [$SUB] Copying Files ... "
      START=$(date -u +%s.%N)


      if test ! -z "${BLOCK}"; then
            # Call the project specific script
            source ${IMPORT_SCRIPT} ${SUB} ${OUTDIR} ${BLOCK}
      else
            # Call the project specific script
            source ${IMPORT_SCRIPT} ${SUB} ${OUTDIR}
      fi

      if ! test -f "${OUTDIR}/rev_phase.nii"; then
            # Failed to import fieldmap, thus ensure this step is not performed
            echo "No Fieldmap has been imported. No geometric distortion will be performed."
            DO_FMAP=0
      fi


      END=$(date -u +%s.%N)
      DIFF=`echo "( $END - $START )" | bc`
      printf "DONE [%.1f s]\n" $DIFF
fi


# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------



# PERFORM Motion correction and other spatial processing in native space
# The order in which those steps are performed is a topic of debate.
# As far as I know, there is no objectively best way of defining this

# Here, I am following the sequence proposed by Jo et al. 2013
# Effective Preprocessing Procedures Virtually Eliminate Distance-Dependent Motion Artifacts in Resting State FMRI
if [ "$DO_REG" -eq "1" ]; then
      printf "[$SUB] Starting Native space preprocessing ... "
      START=$(date -u +%s.%N)


      #=====================================

      {
            START=$(date -u +%s.%N)

            estimate_func_bias()
            {

                  TN=$1
                  OUTDIR=$2
                  ABIN=$3

                  if (( $TN < 10 )); then
                        FI=00$TN
                  elif (( $TN < 100 )); then
                        FI=0$TN
                  else
                        FI=$TN
                  fi

                  {
                  ${ABIN}/N4BiasFieldCorrection -i ${OUTDIR}/tmp/func_data.${FI}.nii -s 4  -o [ ${OUTDIR}/tmp/bcorr_${FI}.nii,${OUTDIR}/tmp/biasfield_${FI}.nii ]
                  } &> /dev/null
            }
            export -f estimate_func_bias

            skull_strip()
            {

                  TN=$1
                  OUTDIR=$2
                  ABIN=$3
                  DO_NLM=$4



                  if (( $TN < 10 )); then
                        FI=00$TN
                  elif (( $TN < 100 )); then
                        FI=0$TN
                  else
                        FI=$TN
                  fi

                  {
                        3dAutomask -apply_prefix ${OUTDIR}/tmp/func_data_n.${FI}.nii -prefix ${OUTDIR}/tmp/mask.${FI}.nii ${OUTDIR}/tmp/func_data.${FI}.nii
                        # In high SNR scans, acquisition noise should approach a Gaussian distribution (find again the REF)
                        # If this is not the case, the the noise is Rician
                        if [ "$DO_NLM" -eq "1" ]; then
                              ${ABIN}/DenoiseImage -x ${OUTDIR}/tmp/mask.${FI}.nii -s 1 -n Rician -i ${OUTDIR}/tmp/func_data_n.${FI}.nii  -o [ ${OUTDIR}/tmp/func_data_nm.${FI}.nii ]
                              mv ${OUTDIR}/tmp/func_data_nm.${FI}.nii ${OUTDIR}/tmp/func_data_n.${FI}.nii
                        fi
                  } &> /dev/null

            }
            export -f skull_strip


            print_debug 'Performing volume-wise automasking'

            NVOLS=`fslnvols ${OUTDIR}/func_data.nii`
            NZ=`3dinfo -nk ${OUTDIR}/func_data.nii`

            print_debug "NVOLS = ${NVOLS}"


            # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            #
            #  Bias field correction
            #
            #  It is highly unlikely that this step is necessary, but in some scanners
            #  the smooth difference in intensities might affect motion correction
            #  I still need to test this, though
            #
            # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            PREF=''

            if test ! -f "${OUTDIR}/slice_acq.txt"; then
                  # File is not present (normally this should happen because data was not imported from DICOM)
                  python ./util/write_slice_timing.py -tr ${TR} -nsl ${NZ} -acq ${ACQ_TYPE} -out ${OUTDIR}/slice_acq.txt
            fi

            # Slomoco is better performed, it seems, prior to anything else
            if [ "$MOCO" == "slomoco" ]; then
                  print_debug 'DOING SLOMOCO'
                  source ./align_pestica.sh ${OUTDIR} ${PREF}func_data ${TR} 1


                  3dresample -prefix ${OUTDIR}/tmp.nii -orient lpi -input ${OUTDIR}/${PREF}func_data.nii
                  3drefit -deoblique ${OUTDIR}/tmp.nii

                  #TODO Remove this line [only here for comparison]
                  3dvolreg -heptic -prefix ${OUTDIR}/r2${PREF}func_data.nii -base 0 -rot_thresh 0.01 -delta 2 \
                        -x_thresh 0.01 -zpad 10 -maxite 60 -1Dfile ${OUTDIR}/motion_estimate.par -maxdisp1D ${OUTDIR}/maximum_disp.1d \
                         ${OUTDIR}/tmp.nii

                  rm -f ${OUTDIR}/tmp.nii


                  3drefit -deoblique ${OUTDIR}/r2${PREF}func_data.nii
                  1dplot -thick -volreg -png ${OUTDIR}/motion_estimate.png -one ${OUTDIR}/motion_estimate.par
                  1dplot -thick -png ${OUTDIR}/maximum_disp.png -one ${OUTDIR}/maximum_disp.1d
                  1dplot -thick -png ${OUTDIR}/maximum_disp_delt.png -one ${OUTDIR}/maximum_disp.1d_delt


                  python ${QCDIR}/save_image_diff.py -o ${OUTDIR} -i ${OUTDIR} -a r2${PREF}func_data.nii -b rp${PREF}func_data.nii -msg1 'PESTICA+3dvorleg' -msg2 'PESTICA+SLOMOCO'
                  cp ${OUTDIR}/slomoco4/*.txt ${OUTDIR}/
                  cp ${OUTDIR}/slomoco4/*.1D ${OUTDIR}/
                  cp ${OUTDIR}/slomoco4/*.png ${OUTDIR}/

                  PREF=rp${PREF}
                  log_command_div 'MOTION CORRECTION'
                  cd ${WDIR}
            fi


            # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            if [ "$DO_SLC" -eq "1" ]; then
                  print_debug 'Slice timing correction'
                  # Slice timing correction  [Parker et al. 2017 -- 10.1016/j.media.2016.08.006]
                  CF=`echo "(1/${TR})/2" | bc -l`

                  if test ! -f "${OUTDIR}/slice_acq.txt"; then
                        # File is not present (normally this should happen because data was not imported from DICOM)
                        python ./util/write_slice_timing.py -tr ${TR} -nsl ${NZ} -acq ${ACQ_TYPE} -out ${OUTDIR}/slice_acq.txt
                  fi


                  if (( $(echo "$TR >= 2.5" | bc -l) )); then
                        # Used to preserve the full power spectrum
                        CF=0.21
                  fi

                  # Correct slice acquisition to the middle of the volume acquisition
                  # To correct to the first slice, set REF_TIME to 0
                  # TODO [10.09.19] Pass this to the parameter section of the file
                  REF_TIME=`echo "${TR}/2" | bc -l`


                  print_debug 'Slice timing correction' "filtershift --in=${OUTDIR}/${PREF}func_data.nii --itl=${ACQ_TYPE} \
                    --cf=${CF} --rt=${REF_TIME} --TR=${TR} --timing=${OUTDIR}/slice_acq.txt --out=${OUTDIR}/a${PREF}func_data.nii"

                  filtershift --in=${OUTDIR}/${PREF}func_data.nii --itl=${ACQ_TYPE}  --cf=${CF} --rt=${REF_TIME} --TR=${TR} --out=${OUTDIR}/a${PREF}func_data.nii
                  gunzip -f ${OUTDIR}/a${PREF}func_data.nii.gz

                  PREF=a${PREF}
                  log_command_div 'Slice Timing correction'
            fi

            if [ "$DO_DEOBL" -eq "1" ]; then
                  print_debug 'Deobliquing volumes'

                  #TODO CHECK AND CORRECT ORIENTATION INFO!

                  #3dWarp -deoblique -newgrid ${DEOBL_VOX} -NN -prefix ${OUTDIR}/w${PREF}func_data.nii ${OUTDIR}/${PREF}func_data.nii
                  #cp ${OUTDIR}/${PREF}func_data.nii ${OUTDIR}/w${PREF}func_data.nii
                  3dresample -prefix ${OUTDIR}/w${PREF}func_data.nii -orient lpi -input ${OUTDIR}/${PREF}func_data.nii
                  3drefit -deoblique ${OUTDIR}/w${PREF}func_data.nii

                  PREF=w${PREF}
                  log_command_div 'DEOBLIQUE'
            fi


            # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


            if [ "${DO_FUNC_BIAS}" -eq "1" ]; then
                  log_command_div 'Functional Bias Correction'

                  mkdir ${OUTDIR}/tmp
                  3dTsplit4D -prefix ${OUTDIR}/tmp/func_data.nii -keep_datum ${OUTDIR}/${PREF}func_data.nii
                  parallel -j4 --line-buffer estimate_func_bias ::: $(seq 0 ${NVOLS}) ::: ${OUTDIR} ::: ${ABIN}

                  3dTcat -prefix ${OUTDIR}/func_biasfield.nii -tr ${TR} ${OUTDIR}/tmp/biasfield_*.nii

                  fslmaths ${OUTDIR}/func_biasfield.nii -Tmean ${OUTDIR}/mean_biasfield.nii
                  fslmaths ${OUTDIR}/${PREF}func_data.nii -div ${OUTDIR}/mean_biasfield.nii ${OUTDIR}/b${PREF}func_data.nii
                  gunzip -f ${OUTDIR}/b${PREF}func_data.nii.gz

                  rm -f ${OUTDIR}/func_biasfield.nii
                  rm -r ${OUTDIR}/tmp
                  #python ${QCDIR}/save_image_diff.py -o ${OUTDIR} -i ${OUTDIR} -a func_data.nii -b ${PREF}func_data.nii -type mean -msg1 'With Bias' -msg2 'Without Bias'
                  log_command_div 'Functional Bias Correction'
                  PREF=b${PREF}
            fi


            # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            # Skull stripping
            # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


            log_command_div 'Skull strpping [3dAutomask on individual volumes]'
            mkdir ${OUTDIR}/tmp
            3dTsplit4D -prefix ${OUTDIR}/tmp/func_data.nii -keep_datum ${OUTDIR}/${PREF}func_data.nii
            parallel -j4 --line-buffer skull_strip ::: $(seq 0 ${NVOLS}) ::: ${OUTDIR} ::: ${ABIN} ::: ${DO_NLM}

            PREF=m${PREF}
            3dTcat -prefix ${OUTDIR}/${PREF}func_data.nii -tr ${TR} ${OUTDIR}/tmp/func_data_n*.nii


            rm -r ${OUTDIR}/tmp
            log_command_div 'Skull stripping'


            # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            # If movement is extreme, this might help. Then again, it is questionable when exactly to do this.
            # At the moment, I'm leaning towards not performing this step and rather censoring offending volumes with the 3dTproject function

            if [ "$DO_DPK" -eq "1" ]; then
                  print_debug 'Despiking [-localedit -NEW]'
                  3dDespike -nomask -NEW -localedit -prefix ${OUTDIR}/d${PREF}func_data.nii ${OUTDIR}/${PREF}func_data.nii
                  python ${QCDIR}/save_image_diff.py -o ${OUTDIR} -i ${OUTDIR} -a ${PREF}func_data.nii -b d${PREF}func_data.nii -msg1 'Before Despike' -msg2 'After Despike'
                  PREF=d${PREF}
                  log_command_div 'Despike'
            fi


            # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            #
            # Motion correction
            #
            # There are 2 methods here: 3dvolreg or slomoco
            # My tests seem to confirm the idea that slomoco performs better than 3dvolreg
            # This procedure, however, is extremely computationaly intensive (i.e. it takes hours instead of minutes to complete)
            #
            # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            if [ "$MOCO" == "3dvolreg" ]; then
                  VOL2MED_COR=`3dTqual -automask ${OUTDIR}/${PREF}func_data.nii`
                  MIN_VAL=100
                  REF_VOL=0
                  IDX=0
                  while IFS=' \t' read -r  qualVal;
                  do
                        if (( $(echo "${qualVal} < ${MIN_VAL}" | bc -l) )); then
                              MIN_VAL=${qualVal}
                              echo 'New MIN : [' $IDX ']: ' ${qualVal}
                              REF_VOL=$IDX

                        fi
                        IDX=$((IDX+1))
                  done <<< "${VOL2MED_COR}";

                  print_debug 'Motion Correction' "3dvolreg -heptic -prefix ${OUTDIR}/r${PREF}func_data.nii -base ${REF_VOL} -rot_thresh 0.02  \
                              -x_thresh 0.01 -zpad 5 -maxite 40 -1Dfile ${OUTDIR}/motion_estimate.par -maxdisp1D ${OUTDIR}/maximum_disp.1d  \
                               ${OUTDIR}/${PREF}func_data.nii"


                  # Motion correction
                  # Performing the realignment to the mean volume increases the processing time quite a lot, without any major benefit
                  # (at least that is what Oakes et al. 2005 say ['Comparison of fMRI motion correction software tools'], though connectivity was not investigated )
                  3dTcat -prefix ${OUTDIR}/ref_vol.nii -TR ${TR} ${OUTDIR}/${PREF}func_data.nii[${REF_VOL}]
                  3dmerge -prefix ${OUTDIR}/ref_vol_smooth.nii -1blur_fwhm 2 ${OUTDIR}/ref_vol.nii


                  3dvolreg -heptic -prefix ${OUTDIR}/r${PREF}func_data.nii -base ${OUTDIR}/ref_vol_smooth.nii -rot_thresh 0.001 -delta 0.15 \
                        -x_thresh 0.001 -zpad 10 -maxite 75 -1Dfile ${OUTDIR}/motion_estimate.par -maxdisp1D ${OUTDIR}/maximum_disp.1d \
                         ${OUTDIR}/${PREF}func_data.nii



                  1dplot -thick -volreg -png ${OUTDIR}/motion_estimate_0.png -one ${OUTDIR}/motion_estimate.par
                  1dplot -thick -png ${OUTDIR}/maximum_disp_0.png -one ${OUTDIR}/maximum_disp.1d
                  1dplot -thick -png ${OUTDIR}/maximum_disp_delt_0.png -one ${OUTDIR}/maximum_disp.1d_delt



                  PREF=r${PREF}

                  log_command_div 'MOTION CORRECTION'

            fi


            if [ "$DO_DPK" -eq "2" ]; then
                  print_debug 'Despiking [-localedit -NEW]'
                  3dDespike -nomask -NEW -localedit -prefix ${OUTDIR}/d${PREF}func_data.nii ${OUTDIR}/${PREF}func_data.nii
                  python ${QCDIR}/save_image_diff.py -o ${OUTDIR} -i ${OUTDIR} -a ${PREF}func_data.nii -b d${PREF}func_data.nii -msg1 'Before Despike' -msg2 'After Despike'
                  PREF=d${PREF}
                  log_command_div 'Despike'
            fi

            # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            # IMPORTANT NOTE
            # Aug 3, 2012  07:08 AM | Erik Beall
            # RE: Is motion correction needed prior PESTICA?
            # Hi Yaroslav,
            # Running the estimation before doing any processing of the data is the default, but I apply the correction _after_ motion correction, because it should have a slightly improved correction (see T Jones et al, Integration of motion correction and physiological noise regression in fMRI. Neuroimage. 2008; 42(2): 582–590. doi: 10.1016/j.neuroimage.2008.05.019), but the effect is minimal.
            # As to running estimation before or after motion correction, I've tried PESTICA estimation using both unprocessed and motion-corrected data for the slicewise ICA stage, with similar results.  I didn't do a detailed group comparison, but I suggest running PESTICA estimation on unprocessed data (before moco) and then apply the resulting estimators with RETROICOR or IRF-RETROICOR correction to volumetric motion-corrected data.

            # !!!!
            # Feb 21, 2013  12:02 PM | Erik Beall
            # RE: Is motion correction needed prior PESTICA?
            # Yaroslav, I don't know how you ended up integrating PESTICA into your pipeline, but I've put out a new release and I wanted to point out that the new version no longer allows you to apply PESTICA to a different datafile than what you use as input.  In practice, you can run PESTICA either before or after motion correction, but I recommend after.  I've learned that applying PESTICA estimators from a differently-corrected dataset onto another corrected dataset of the same raw data sometimes gives really bad results.  Doing it the default way (in v2.0) seems to work very robustly in nealry all circumstances.
            # Erik

            if [ "$DO_PEST" -eq "1" ] && [ ! "$MOCO" == "slomoco" ]; then
                  print_debug 'PESTICA4 cardiac/respiratory effects correction'
                  # 0 here indicates that SLOMOCO procedure should not be done
                  source ./align_pestica.sh ${OUTDIR} ${PREF}func_data ${TR} 0
                  python ${QCDIR}/save_image_diff.py -o ${OUTDIR} -i ${OUTDIR} -a ${PREF}func_data.nii -b p${PREF}func_data.nii -msg1 'Before PESTICA' -msg2 'After PESTICA'
                  PREF=p${PREF}
                  log_command_div 'Physiological Signal Correction'
            fi



            # Create native space brain mask
            3dAutomask -prefix ${OUTDIR}/nat_mask.nii ${OUTDIR}/${PREF}func_data.nii



            # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


            # Distortion correction using reverse phase acquisition
            if [ "$DO_FMAP" -eq "1" ]; then
                  print_debug 'Reverse phase fieldmap estimation and distortion correction'

                  # NOTE: Topup fails with odd number of slices, thus we zero pad if this occurs
                  #https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=fsl;67dcb45c.1209
                  X=`3dinfo -ni ${OUTDIR}/${PREF}func_data.nii`
                  Y=`3dinfo -nj ${OUTDIR}/${PREF}func_data.nii`
                  Z=`3dinfo -nk ${OUTDIR}/${PREF}func_data.nii`


                  HAS_ODD=0
                  if (( $Z % 2 )); then
                        # Z is odd
                        Z=$((Z+1))
                        HAS_ODD=1
                  fi

                  if (( $X % 2 )); then
                        # Z is odd
                        X=$((X+1))
                        HAS_ODD=1
                  fi

                  if (( $Y % 2 )); then
                        # Z is odd
                        Y=$((Y+1))
                        HAS_ODD=1
                  fi

                  if [ "${HAS_ODD}" -eq "1" ]; then
                        fslroi ${OUTDIR}/${PREF}func_data ${OUTDIR}/tmp 0 $X 0 $Y 0 $Z
                        gunzip -f ${OUTDIR}/tmp.nii.gz
                        mv ${OUTDIR}/tmp.nii  ${OUTDIR}/${PREF}func_data.nii

                        fslroi ${OUTDIR}/rev_phase ${OUTDIR}/tmp 0 $X 0 $Y 0 $Z
                        gunzip -f ${OUTDIR}/tmp.nii.gz
                        mv ${OUTDIR}/tmp.nii  ${OUTDIR}/rev_phase.nii
                  fi

                  # Use 5 volumes with normal and reverse phase encoding and create brain mask
                  HVOL=$((NVOLS/2))
                  VI=$((HVOL-2))
                  VF=$((HVOL+2))
                  3dTcat -tr ${TR} -prefix ${OUTDIR}/blipdown.nii ${OUTDIR}/${PREF}func_data.nii[${VI}..${VF}]
                  3dTcat -tr ${TR} -prefix ${OUTDIR}/blipup.nii ${OUTDIR}/rev_phase.nii

                  3dAutomask -apply_prefix ${OUTDIR}/blipup_m.nii ${OUTDIR}/blipup.nii
                  3dAutomask -apply_prefix ${OUTDIR}/blipdown_m.nii ${OUTDIR}/blipdown.nii

                  # Concatenate both PE to create single 10-volume series from where the fieldmap will be esimated
                  3dTcat -tr ${TR} -prefix ${OUTDIR}/topup_data.nii ${OUTDIR}/blipdown_m.nii ${OUTDIR}/blipup_m.nii


                  rm -f ${OUTDIR}/blipdown.nii
                  rm -f ${OUTDIR}/blipup.nii
                  rm -f ${OUTDIR}/blipdown_m.nii
                  rm -f ${OUTDIR}/blipup_m.nii


                  python ./util/write_topup_file.py -out ${OUTDIR}/topupfield.txt

                  # Fieldmap estimation [ONLY estimation]
                  print_debug 'TOPUP command' "topup --imain=${OUTDIR}/topup_data.nii --datain=${OUTDIR}/topupfield.txt --out=${OUTDIR}/topup_results --fout=${OUTDIR}/fieldmap --iout=${OUTDIR}/w_topup_data --estmov=0,0,0  --regmod=membrane_energy  --minmet=1,1,1 --verbose --warpres=8,6,4 --miter=32,12,4 --subsamp=4,2,1 --fwhm=6,4,2"

                  topup --imain=${OUTDIR}/topup_data.nii --datain=${OUTDIR}/topupfield.txt --out=${OUTDIR}/topup_results --fout=${OUTDIR}/fieldmap --iout=${OUTDIR}/w_topup_data \
                        --estmov=0,0,0  --regmod=membrane_energy  --minmet=1,1,1 --verbose --warpres=8,6,4 --miter=32,12,4 --subsamp=2,2,1 --fwhm=6,4,2

                  applytopup --imain=${OUTDIR}/${PREF}func_data.nii --datain=${OUTDIR}/topupfield.txt --topup=${OUTDIR}/topup_results --inindex=1 --method=jac --interp=spline --out=${OUTDIR}/u${PREF}func_data
                  gunzip -f ${OUTDIR}/u${PREF}func_data.nii.gz
                  3dAutomask -apply_prefix ${OUTDIR}/tmp.nii ${OUTDIR}/u${PREF}func_data.nii
                  rm ${OUTDIR}/u${PREF}func_data.nii
                  mv ${OUTDIR}/tmp.nii ${OUTDIR}/u${PREF}func_data.nii
                  python ${QCDIR}/save_image_diff.py -o ${OUTDIR} -i ${OUTDIR} -a ${PREF}func_data.nii -b u${PREF}func_data.nii -type mean -msg1 'Before FMAP' -msg2 'After FMAP'
                  PREF=u${PREF}
                  log_command_div 'FIELDMAP CORRECTION'

            fi


            # Distortion correction using acquired fieldmap (magnitude and phase)
            if [ "$DO_FMAP" -eq "2" ]; then
                  print_debug 'Fieldmap correction'

                  fslmaths ${OUTDIR}/${PREF}func_data -Tmean  ${OUTDIR}/mean_func


                  flirt  -in ${OUTDIR}/fmap_mag_brain -ref ${OUTDIR}/mean_func -omat ${OUTDIR}/fmap2func

                  flirt  -in ${OUTDIR}/fmap_mag -ref ${OUTDIR}/mean_func -applyxfm -init ${OUTDIR}/fmap2func -out ${OUTDIR}/fmap_mag_func
                  flirt  -in ${OUTDIR}/fmap_phase_rads -ref ${OUTDIR}/mean_func -applyxfm -init ${OUTDIR}/fmap2func -out ${OUTDIR}/fmap_phase_rads_func


                  fugue --loadfmap=${OUTDIR}/fmap_phase_rads_func -s 1 --despike --savefmap=${OUTDIR}/fmap_phase_rads_func

                  #TODO Pass unwarpdir to configuration
                  fugue -i ${OUTDIR}/${PREF}func_data.nii --dwell=${EES} --loadfmap=${OUTDIR}/fmap_phase_rads_func -u ${OUTDIR}/u${PREF}func_data.nii --unwarpdir=y-

                  gunzip -f ${OUTDIR}/u${PREF}func_data.nii.gz
                  python ${QCDIR}/save_image_diff.py -o ${OUTDIR} -i ${OUTDIR} -a ${PREF}func_data.nii -b u${PREF}func_data.nii -type mean -msg1 'Before FMAP' -msg2 'After FMAP'

                  PREF=u${PREF}
                  log_command_div 'FIELDMAP CORRECTION'
            fi


            # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



            if [ "$DO_ING" -eq "1" ]; then
                  print_debug 'Global intensity normalization' 'fslmaths ${OUTDIR}/${PREF}func_data -ing ${GLOB_VAL} ${OUTDIR}/i${PREF}func_data'
                  fslmaths ${OUTDIR}/${PREF}func_data -ing ${GLOB_VAL} ${OUTDIR}/i${PREF}func_data
                  gunzip -f ${OUTDIR}/i${PREF}func_data.nii.gz
                  PREF=i${PREF}
            fi


            # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            # Update brain mask
            rm ${OUTDIR}/nat_mask.nii
            3dAutomask -prefix ${OUTDIR}/nat_mask.nii -apply_prefix ${OUTDIR}/proc_data_native.nii ${OUTDIR}/${PREF}func_data.nii

            if [ "$DO_CLEAN" -eq "1" ]; then
                  cp ${OUTDIR}/pestica4/*.png ${OUTDIR}/
                  rm -r -f ${OUTDIR}/pestica4

                  cp ${OUTDIR}/slomoco4/*.png ${OUTDIR}/
                  cp ${OUTDIR}/slomoco4/*.1D ${OUTDIR}/
                  cp ${OUTDIR}/slomoco4/*.txt ${OUTDIR}/
                  rm -r -f ${OUTDIR}/slomoco4



                  rm -f ${OUTDIR}/*func_data.nii
                  rm -f ${OUTDIR}/*bias*.nii
                  rm -f ${OUTDIR}/*func_data.nii.gz
                  rm -f ${OUTDIR}/*.BRIK
                  rm -f ${OUTDIR}/*.HEAD
            fi


      } &> ${OUTDIR}/01_moco.log
      END=$(date -u +%s.%N)
      DIFF=`echo "( $END - $START )" | bc`
      printf "DONE [%.1f s]\n" $DIFF

fi




# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# END OF NATIVE VOLUMETRIC PROCESSING
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------


# TODO Add the possibility to import this from somewhere
# This would be useful when there are multiple functional blocks (no need process anatomical and register to MNI again in those cases)

if [ "${EXTRACT_SURF}" -eq "1" ]; then
    matlab "-nodesktop -nosplash " <<<"cd '${WDIR}'; anat_proc_surf('${OUTDIR}', 't1_f.nii'); exit;"
fi

# Calculate WARP to group template
if [ "$DO_REG2" -eq "1" ]; then
      printf "[$SUB] Functional <-> Anat <-> MNI Registration ... "
      START=$(date -u +%s.%N)


      {

      # Bias field correction and Gaussian noise reduction
      # TODO Pass to parameter section
      # Necessary files must be copied in case this section does not run (e.g. multiple blocks of functional data for the same participant)
      DO_ANAT=1


      # Note: slomoco will run 3dAllineate, which removes oblique information from the header
      # As such, it is a good idea to deoblique the T1 volume in order to match orientation
      if [ "$DO_DEOBL" -eq "1" ]; then
            #3drefit -deoblique ${OUTDIR}/t1_f.nii
            3dresample -prefix ${OUTDIR}/t1_f.nii -orient lpi -input ${OUTDIR}/t1.nii

      else
            cp ${OUTDIR}/t1.nii ${OUTDIR}/t1_f.nii
      fi



      ${ABIN}/N4BiasFieldCorrection -i ${OUTDIR}/t1_f.nii -s 4  -o [ ${OUTDIR}/t1_fb.nii,${OUTDIR}/biasfield.nii ]


      if [ "$BREXT_TYPE" -eq "1" ]; then
            # [NOTE]
            # deepbrain (used here) requires an older version of tensorflow (v1x)
            # the way araound this is to edit lines 18, 19, and 20 (all with tf) in file
            # /home/luna.kuleuven.be/u0101486/Software/anaconda3/lib/python3.7/site-packages/deepbrain/extractor.py
            # "According to TF 1:1 Symbols Map, in TF 2.0 you should use tf.compat.v1.Session() instead of tf.Session()"
            # Source: https://stackoverflow.com/questions/55142951/tensorflow-2-0-attributeerror-module-tensorflow-has-no-attribute-session
            python ${WDIR}/extract_brain.py -i "${OUTDIR}/t1_fb.nii" -o "${OUTDIR}/AnatMask.nii"
      fi


      if [ "$BREXT_TYPE" -eq "2" ]; then
            #TODO Add FSL's bet extraction procedure
            echo 'TODO'
      fi

      3dcalc -a ${OUTDIR}/t1_fb.nii -b ${OUTDIR}/AnatMask.nii -expr 'a*b' -prefix ${OUTDIR}/t1_brain.nii
      ${ABIN}/DenoiseImage -i ${OUTDIR}/t1_brain.nii  -o [ ${OUTDIR}/t1_frn.nii,${OUTDIR}/t1_noise.nii ]

      # [NOTE]
      #  BBR procedure likely fails when tissue contrast is poor
      3dUnifize -prefix ${OUTDIR}/anat_ref.nii ${OUTDIR}/t1_frn.nii

      # Anat2MNI. The results of this transform are also used to obtain priors for segmentation
      SRC=${OUTDIR}/anat_ref.nii

      ${ABIN}/antsRegistration -d 3 -r [$MNI_REF,$SRC,0] -v 1 \
                  -m MI[$MNI_REF,$SRC,1,64] -t translation[0.2] -c [500x50,1.e-8,20] \
                  -s 6x4vox -f 3x2 -l 1 -n BSpline \
                  -m MI[$MNI_REF,$SRC,1,64,Regular,0.5] -t rigid[0.2] -c [500x50,1.e-8,20] \
                  -s 6x4vox -f 3x2 -l 1 -n BSpline \
                  -m MI[$MNI_REF,$SRC,1,64,Regular,0.5] -t affine[0.2] -c [500x200x200,1.e-8,10] \
                  -s 6x6x4vox -f 4x3x2 -l 1 -n BSpline \
                  -m CC[$MNI_REF,$SRC,1,3] -t SyN[0.2,3] -c [220x150x25,1.e-7,10] \
                  -s 6x3x0vox -f 4x2x1 -l 1 -n BSpline \
                  -o [${OUTDIR}/anat2group,${OUTDIR}/anat2group.nii]


      if [ "$USE_TISSUE_PRIOR" -eq "1" ]; then

            ${ABIN}/antsApplyTransforms \
                  -i ${GM_PRIOR} \
                  --float \
                  -r ${SRC} \
                  -o ${OUTDIR}/gm_tpm_anat.nii \
                  -t [${OUTDIR}/anat2group0GenericAffine.mat,1] \
                  -t [${OUTDIR}/anat2group1InverseWarp.nii.gz,0] \
                  -n Linear

            ${ABIN}/antsApplyTransforms \
                  -i ${WM_PRIOR} \
                  --float \
                  -r ${SRC} \
                  -o ${OUTDIR}/wm_tpm_anat.nii \
                  -t [${OUTDIR}/anat2group0GenericAffine.mat,1] \
                  -t [${OUTDIR}/anat2group1InverseWarp.nii.gz,0] \
                  -n Linear

            ${ABIN}/antsApplyTransforms \
                  -i ${CSF_PRIOR} \
                  --float \
                  -r ${SRC} \
                  -o ${OUTDIR}/csf_tpm_anat.nii \
                  -t [${OUTDIR}/anat2group0GenericAffine.mat,1] \
                  -t [${OUTDIR}/anat2group1InverseWarp.nii.gz,0] \
                  -n Linear
            fast --nobias --segments --verbose -A ${OUTDIR}/gm_tpm_anat.nii ${OUTDIR}/wm_tpm_anat.nii ${OUTDIR}/csf_tpm_anat.nii ${OUTDIR}/t1_frn.nii
      else
            fast --nobias --segments --verbose ${OUTDIR}/t1_frn.nii
      fi

      gunzip -f ${OUTDIR}/t1_frn_seg_0.nii.gz
      mv ${OUTDIR}/t1_frn_seg_0.nii ${OUTDIR}/csf_mask.nii

      gunzip -f ${OUTDIR}/t1_frn_seg_1.nii.gz
      mv ${OUTDIR}/t1_frn_seg_1.nii ${OUTDIR}/gm_mask.nii

      gunzip -f ${OUTDIR}/t1_frn_seg_2.nii.gz
      mv ${OUTDIR}/t1_frn_seg_2.nii ${OUTDIR}/wm_mask.nii

      ${ABIN}/antsApplyTransforms \
            -i ${OUTDIR}/gm_mask.nii \
            --float \
            -r ${MNI_REF_2mm} \
            -o ${OUTDIR}/gm_mask_mni.nii \
            -t [${OUTDIR}/anat2group1Warp.nii.gz,0] \
            -t [${OUTDIR}/anat2group0GenericAffine.mat,0] \
            -n Linear

      ${ABIN}/antsApplyTransforms \
            -i ${OUTDIR}/wm_mask.nii \
            --float \
            -r ${MNI_REF_2mm} \
            -o ${OUTDIR}/wm_mask_mni.nii \
            -t [${OUTDIR}/anat2group1Warp.nii.gz,0] \
            -t [${OUTDIR}/anat2group0GenericAffine.mat,0] \
            -n Linear

      ${ABIN}/antsApplyTransforms \
            -i ${OUTDIR}/csf_mask.nii \
            --float \
            -r ${MNI_REF_2mm} \
            -o ${OUTDIR}/csf_mask_mni.nii \
            -t [${OUTDIR}/anat2group1Warp.nii.gz,0] \
            -t [${OUTDIR}/anat2group0GenericAffine.mat,0] \
            -n Linear

      # Prepare reference functional image
      # Strictly speaking, it is not necessary to denoise, but my tests indicate that this might help in some cases
      # I cannot see a reason it would negatively impact anything
      fslmaths "${OUTDIR}/proc_data_native.nii" -Tmedian -thr 0 "${OUTDIR}/ref_func.nii"
      gunzip -f ${OUTDIR}/ref_func.nii.gz
      ${ABIN}/N4BiasFieldCorrection -i ${OUTDIR}/ref_func.nii -s 1  -o [ ${OUTDIR}/ref_func_bf.nii ]
      #${ABIN}/DenoiseImage -i ${OUTDIR}/ref_func_b.nii  -o [ ${OUTDIR}/ref_func_bf.nii]
      rm -f ${OUTDIR}/ref_func_bfm.nii
      rm -f ${OUTDIR}/ref_func_bfmu.nii
      3dAutomask  -apply_prefix ${OUTDIR}/ref_func_bfm.nii ${OUTDIR}/ref_func_bf.nii
      3dUnifize -prefix ${OUTDIR}/ref_func_bfmu.nii -EPI -GM ${OUTDIR}/ref_func_bfm.nii

      3dresample -prefix ${OUTDIR}/ref_func_bfmu_3.nii -dxyz 3 3 3 -rmode 'Cubic' -input ${OUTDIR}/ref_func_bfmu.nii
      ${ABIN}/DenoiseImage -i ${OUTDIR}/ref_func_bfmu_3.nii  -o [ ${OUTDIR}/ref_func_bfmu_3d.nii]

      3dresample -prefix ${OUTDIR}/ref_func_bfmu_2.nii -dxyz 2 2 2 -rmode 'Cubic' -input ${OUTDIR}/ref_func_bfmu_3d.nii
      ${ABIN}/DenoiseImage -i ${OUTDIR}/ref_func_bfmu_2.nii  -o [ ${OUTDIR}/ref_func_bfmu_2d.nii]

      mv ${OUTDIR}/ref_func_bfmu_2d.nii ${OUTDIR}/ref_func_bfm.nii


      if [ "${F2A_FUNC}" = "flirt_bbr" ]; then
        flirt -in ${OUTDIR}/ref_func_bfm.nii -out ${OUTDIR}/afunc2anat.nii -ref ${OUTDIR}/anat_ref.nii \
            -cost bbr -omat ${OUTDIR}/afunc2anat.mat -dof 6 \
            -wmseg ${OUTDIR}/wm_mask.nii -bbrtype local_abs
        gunzip -f ${OUTDIR}/afunc2anat.nii.gz
        ${C3DPATH}/c3d_affine_tool -ref ${OUTDIR}/t1_frn.nii -src ${OUTDIR}/ref_func_bf.nii ${OUTDIR}/afunc2anat.mat -fsl2ras -oitk ${OUTDIR}/func2anat.mat
      fi

      if [ "${F2A_FUNC}" = "flirt" ]; then
          flirt -in ${OUTDIR}/ref_func_bfm.nii -out ${OUTDIR}/afunc2anat.nii -ref ${OUTDIR}/anat_ref.nii \
              -cost normmi -omat ${OUTDIR}/afunc2anat.mat  \
              -dof 6
          gunzip -f ${OUTDIR}/afunc2anat.nii.gz
          ${C3DPATH}/c3d_affine_tool -ref ${OUTDIR}/t1_frn.nii -src ${OUTDIR}/ref_func_bf.nii ${OUTDIR}/afunc2anat.mat -fsl2ras -oitk ${OUTDIR}/func2anat.mat
      fi

      if [ "${F2A_FUNC}" = "afni_twopass" ]; then
        3dAllineate -input ${OUTDIR}/ref_func_bfm.nii -base ${OUTDIR}/anat_ref.nii \
            -prefix ${OUTDIR}/afunc2anat.nii -1Dmatrix_save ${OUTDIR}/func2anat_afni.mat \
            -cmass -cost lpc+hel -warp shift_rotate -interp cubic -nmatch '30%'
            python ${WDIR}/util/afni2ras_affine.py -in ${OUTDIR}/func2anat_afni.aff12.1D -out ${OUTDIR}/afunc2anat.mat
            rm ${OUTDIR}/func2anat_afni.aff12.1D
            ${C3DPATH}/c3d_affine_tool -ref ${OUTDIR}/t1_frn.nii -src ${OUTDIR}/ref_func_bf.nii ${OUTDIR}/afunc2anat.mat -oitk ${OUTDIR}/func2anat.mat
      fi

      if [ "${F2A_FUNC}" = "afni_onepass" ]; then
        rm -f ${OUTDIR}/afunc2anat.nii
        3dAllineate -input ${OUTDIR}/ref_func_bfm.nii -base ${OUTDIR}/anat_ref.nii \
            -prefix ${OUTDIR}/afunc2anat.nii -1Dmatrix_save ${OUTDIR}/func2anat_afni \
            -onepass -cmass -cost lpc+hel -warp shift_rotate -interp cubic -nmatch '30%'

        #cat_matvec ${OUTDIR}/func2anat_afni.aff12.1D -4x4 > ${OUTDIR}/afunc2anat.mat
        python ${WDIR}/util/afni2ras_affine.py -in ${OUTDIR}/func2anat_afni.aff12.1D -out ${OUTDIR}/afunc2anat.mat
        rm ${OUTDIR}/func2anat_afni.aff12.1D
        ${C3DPATH}/c3d_affine_tool -ref ${OUTDIR}/t1_frn.nii -src ${OUTDIR}/ref_func_bf.nii ${OUTDIR}/afunc2anat.mat -oitk ${OUTDIR}/func2anat.mat
      fi

      python ${QCDIR}/QC_check_reg.py -out ${OUTDIR} -in ${OUTDIR}   \
                        -im1 ${OUTDIR}/afunc2anat.nii -im2 ${OUTDIR}/anat_ref.nii -outf 'func2anat' \
                        -dpi 300


      # Try to use anatomical T1 to correct field inhomogeneity.
      # Requires good tissue constrast in the functional data
      # IF using ICA-FIX to clean the data, it is probably a good idea to already apply the transform here.
      # Otherwise, during the nromalisation step is better, as it interpolates only once
      if [ "$DO_ANAT_FMAP" -eq "1" ]; then
            SRC=${OUTDIR}/afunc2anat.nii
            REF=${OUTDIR}/t1_frn.nii

            ${ABIN}/antsRegistration -d 3 -v 0 \
                  -m CC[$REF,$SRC,1,3] -t SyN[0.1,3] -c [15,1.e-6,10] -g 0x1x0.1  \
                  -s 3vox -f 2 -l 1 -n BSpline \
                  -o [${OUTDIR}/func2anat,${OUTDIR}/awfunc2anat.nii]

            python ${QCDIR}/QC_check_reg.py -out ${OUTDIR} -in ${OUTDIR}   \
                              -im1 ${OUTDIR}/awfunc2anat.nii -im2 ${OUTDIR}/anat_ref.nii -outf 'func2anat_und' \
                              -dpi 300


            if [ "${CORRECT_ANAT_FMAP_NATIVE}" -eq "1" ]; then
                  apply_anat_fmap()
                  {
                        TN=$1
                        OUTDIR=$2
                        ABIN=$3
                        USBPATH=$4
                        REF=${5}


                        if (( $TN < 10 )); then
                              FI=00$TN
                        elif (( $TN < 100 )); then
                              FI=0$TN
                        else
                              FI=$TN
                        fi

                        ${ABIN}/antsApplyTransforms -v 0 \
                              -i ${OUTDIR}/tmp/proc_data.$FI.nii \
                              --float \
                              -r ${REF} \
                              -o ${OUTDIR}/tmp/proc_data_corr_$FI.nii \
                              -t [${OUTDIR}/func2anat.mat,1] \
                              -t [${OUTDIR}/func2anat0Warp.nii.gz,0] \
                              -t [${OUTDIR}/func2anat.mat,0] \
                              -n BSpline


                        # Remove negative values due to BSpline interpolation (mostly in the background of the images)
                        3dmerge -prefix ${OUTDIR}/tmp/proc_data_corr_thr_$FI.nii -1noneg ${OUTDIR}/tmp/proc_data.$FI.nii

                  }
                  export -f apply_anat_fmap

                  mkdir ${OUTDIR}/tmp
                  NVOLS=`fslnvols ${OUTDIR}/proc_data_native.nii`
                  3dTsplit4D -prefix ${OUTDIR}/tmp/proc_data.nii -keep_datum ${OUTDIR}/proc_data_native.nii
                  parallel -j5 --line-buffer apply_anat_fmap ::: $(seq 0 ${NVOLS}) ::: ${OUTDIR} ::: ${ABIN} ::: ${USBPATH} ::: ${OUTDIR}/ref_func_bfm.nii


                  3dTcat -prefix ${OUTDIR}/proc_data_native_fmap.nii -tr ${TR} ${OUTDIR}/tmp/proc_data_corr_thr_*.nii
                  rm -r ${OUTDIR}/tmp


                  # This here is useful for quality control
                  ${ABIN}/antsApplyTransforms -v 0 \
                        -i ${OUTDIR}/t1_frn.nii \
                        --float \
                        -r ${OUTDIR}/ref_func_bfm.nii \
                        -o ${OUTDIR}/t1_func.nii \
                        -t [${OUTDIR}/func2anat.mat,1] \
                        -n Linear



                  python ~/Documents/rs_proc/QC_funcs/save_fmap_comp.py -o ${OUTDIR} -i ${OUTDIR} -a proc_data_native.nii -b proc_data_native_fmap.nii -c t1_func.nii -msg1 'Before FMAP' -msg2 'After FMAP'
                  rm -f ${OUTDIR}/proc_data_native.nii
                  mv ${OUTDIR}/proc_data_native_fmap.nii ${OUTDIR}/proc_data_native.nii

            fi
      fi




      if [ "$DO_CLEAN" -eq "1" ]; then
            #rm -f ${OUTDIR}/t1.nii
            rm -f ${OUTDIR}/t1_func.nii
            rm -f ${OUTDIR}/t1_f.nii
            rm -f ${OUTDIR}/t1_fb.nii
            rm -f ${OUTDIR}/biasfield.nii
            rm -f ${OUTDIR}/t1_noise.nii
            rm -f ${OUTDIR}/AnatMask.nii
            rm -f ${OUTDIR}/t1_frn_pve_*.nii.gz
            rm -f ${OUTDIR}/t1_frn_pveseg.nii.gz

      fi
      } &> ${OUTDIR}/02_ant_norm.log

      END=$(date -u +%s.%N)
      DIFF=`echo "( $END - $START )" | bc`
      printf "DONE [%.1f s]\n" $DIFF
fi



# Project the functional data onto the subject's estimated brain surface
if [ "$DO_FUNC_SURF" -eq "1" ]; then
      func_to_anat()
      {

            TN=$1
            OUTDIR=$2
            ABIN=$3
            CORRECT_ANAT_FMAP_NATIVE=$4

            if (( $TN < 10 )); then
                  FI=00$TN
            elif (( $TN < 100 )); then
                  FI=0$TN
            else
                  FI=$TN
            fi

            {
            REF=${OUTDIR}/t1_frn.nii
            if [ "${CORRECT_ANAT_FMAP_NATIVE}" -eq "0" ] && test -f "${OUTDIR}/func2anat0Warp.nii.gz"; then
                  ${ABIN}/antsApplyTransforms -v 0 \
                        -i ${OUTDIR}/tmp/proc_data.$FI.nii \
                        --float \
                        -r ${REF} \
                        -o ${OUTDIR}/tmp/proc_data_ana_$FI.nii \
                        -t [${OUTDIR}/func2anat0Warp.nii.gz,0] \
                        -t [${OUTDIR}/func2anat.mat,0] \
                        -n BSpline
            else
                  ${ABIN}/antsApplyTransforms -v 0 \
                        -i ${OUTDIR}/tmp/proc_data.$FI.nii \
                        --float \
                        -r ${REF} \
                        -o ${OUTDIR}/tmp/proc_data_ana_$FI.nii \
                        -t [${OUTDIR}/func2anat.mat,0] \
                        -n BSpline
            fi

            3dresample -input ${OUTDIR}/tmp/proc_data_ana_$FI.nii -prefix ${OUTDIR}/tmp/proc_anaf_MNI_$FI.nii -master ${OUTDIR}/tmp/proc_data.$FI.nii
            } &> /dev/null

      }
      export -f func_to_anat

      mkdir ${OUTDIR}/tmp
      3dTsplit4D -prefix ${OUTDIR}/tmp/func_data.nii -keep_datum ${OUTDIR}/func_data.nii
      parallel -j4 --line-buffer estimate_func_bias ::: $(seq 0 ${NVOLS}) ::: ${OUTDIR} ::: ${ABIN} ::: ${CORRECT_ANAT_FMAP_NATIVE}


      3dTcat -prefix ${OUTDIR}/proc_data_anatf.nii -tr ${TR} ${OUTDIR}/tmp/proc_anaf_MNI_*.nii
      rm -f -r ${OUTDIR}/tmp


fi
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------


if [ "$DO_ATLAS_NAT" -eq "1" ]; then

      printf "[$SUB] Registering Atlases to Native Space ... "
      START=$(date -u +%s.%N)
      {
      #TODO [11.10.19] Change REF location to be configurable
      # This may be unnecessary, though
      REF=${WDIR}/atlas/MNI152_T1_1mm_brain.nii
      AAL=${WDIR}/atlas//aal/aal2/AAL2_NoCRB.nii
      AAL2=${WDIR}/atlas//aal/aal2/AAL2.nii
      LOCAL_GLOBAL=${WDIR}/atlas//local_global/Schaefer2018_400Parcels_17Networks_order_FSLMNI152_1mm.nii
      SRC=${OUTDIR}/t1_frn.nii

      #TODO Re-enable this
      ${ABIN}/antsRegistration -d 3 -r [$REF,$SRC,0] -v 1 \
                  -m MI[$MNI_REF,$SRC,1,32] -t translation[0.1] -c [500,5.e-7,20] \
                  -s 3vox -f 3 -l 1 -n BSpline \
                  -m MI[$MNI_REF,$SRC,1,32,Regular,0.25] -t rigid[0.1] -c [500,5.e-7,20] \
                  -s 3vox -f 3 -l 1 -n BSpline \
                  -m MI[$MNI_REF,$SRC,1,32,Regular,0.25] -t affine[0.1] -c [500x100x10,5.e-7,10] \
                  -s 2x1x0vox -f 3x2x1 -l 1 -n BSpline \
                  -m CC[$MNI_REF,$SRC,1,3] -t SyN[0.1,3] -c [200x50x10,1.e-7,10] \
                  -s 6x2x0vox -f 4x2x1 -l 1 -n BSpline \
                  -o [${OUTDIR}/anat2atlas,${OUTDIR}/anat2atlas.nii]


      if [ "${CORRECT_ANAT_FMAP_NATIVE}" -eq "0" ] && test -f "${OUTDIR}/func2anat0Warp.nii.gz"; then
            ${ABIN}/antsApplyTransforms -v 1 \
                  -i ${AAL} \
                  --float \
                  -r ${OUTDIR}/ref_func_bf.nii \
                  -o ${OUTDIR}/aal2_nat_atlas.nii \
                  -t [${OUTDIR}/func2anat.mat,1] \
                  -t [${OUTDIR}/func2anat0InverseWarp.nii.gz,0] \
                  -t [${OUTDIR}/anat2atlas0GenericAffine.mat,1] \
                  -t [${OUTDIR}/anat2atlas1InverseWarp.nii.gz,0] \
                  -n NearestNeighbor


                  ${ABIN}/antsApplyTransforms -v 1 \
                        -i ${AAL2} \
                        --float \
                        -r ${OUTDIR}/ref_func_bf.nii \
                        -o ${OUTDIR}/aal2_nat_atlas.nii \
                        -t [${OUTDIR}/func2anat.mat,1] \
                        -t [${OUTDIR}/func2anat0InverseWarp.nii.gz,0] \
                        -t [${OUTDIR}/anat2atlas0GenericAffine.mat,1] \
                        -t [${OUTDIR}/anat2atlas1InverseWarp.nii.gz,0] \
                        -n NearestNeighbor


            ${ABIN}/antsApplyTransforms -v 1 \
                  -i ${LOCAL_GLOBAL} \
                  --float \
                  -r ${OUTDIR}/ref_func_bf.nii \
                  -o ${OUTDIR}/local_global400_nat_atlas.nii \
                  -t [${OUTDIR}/func2anat.mat,1] \
                  -t [${OUTDIR}/func2anat0InverseWarp.nii.gz,0] \
                  -t [${OUTDIR}/anat2atlas0GenericAffine.mat,1] \
                  -t [${OUTDIR}/anat2atlas1InverseWarp.nii.gz,0] \
                  -n NearestNeighbor
      else
            ${ABIN}/antsApplyTransforms -v 1 \
                  -i ${LOCAL_GLOBAL} \
                  --float \
                  -r ${OUTDIR}/ref_func_bf.nii \
                  -o ${OUTDIR}/local_global400_nat_atlas.nii \
                  -t [${OUTDIR}/func2anat.mat,1] \
                  -t [${OUTDIR}/anat2atlas0GenericAffine.mat,1] \
                  -t [${OUTDIR}/anat2atlas1InverseWarp.nii.gz,0] \
                  -n NearestNeighbor

            ${ABIN}/antsApplyTransforms -v 1 \
                  -i ${AAL} \
                  --float \
                  -r ${OUTDIR}/ref_func_bf.nii \
                  -o ${OUTDIR}/aal2_nat_atlas.nii \
                  -t [${OUTDIR}/func2anat.mat,1] \
                  -t [${OUTDIR}/anat2atlas0GenericAffine.mat,1] \
                  -t [${OUTDIR}/anat2atlas1InverseWarp.nii.gz,0] \
                  -n NearestNeighbor

            ${ABIN}/antsApplyTransforms -v 1 \
                  -i ${AAL2} \
                  --float \
                  -r ${OUTDIR}/ref_func_bf.nii \
                  -o ${OUTDIR}/aal2_nat_atlas.nii \
                  -t [${OUTDIR}/func2anat.mat,1] \
                  -t [${OUTDIR}/anat2atlas0GenericAffine.mat,1] \
                  -t [${OUTDIR}/anat2atlas1InverseWarp.nii.gz,0] \
                  -n NearestNeighbor
      fi
      } &> ${OUTDIR}/02b_atlas_norm.log

      END=$(date -u +%s.%N)
      DIFF=`echo "( $END - $START )" | bc`
      printf "DONE [%.1f s]\n" $DIFF

fi



# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------


if  [ "$EXTRACT_NUIS" -eq "1" ]; then

      printf "[$SUB] Extracting Nuisance Regressors ... "
      START=$(date -u +%s.%N)
      {

      echo "Using ${REG_MODEL} model for nuisance regression"

      # Warp tissue maps from anatomical to functional
      # If anatomical 'FMAP' correction was not done, only apply affine transform
      if [ "${DO_ANAT_FMAP}" == "1" ]; then

            ${ABIN}/antsApplyTransforms -v 1 \
                  -i ${OUTDIR}/csf_mask.nii \
                  --float \
                  -r ${OUTDIR}/ref_func_bf.nii \
                  -o ${OUTDIR}/csf_mask_nat.nii \
                  -t [${OUTDIR}/func2anat.mat,1] \
                  -t [${OUTDIR}/func2anat0InverseWarp.nii.gz,0] \
                  -n Linear

            ${ABIN}/antsApplyTransforms -v 1 \
                  -i ${OUTDIR}/wm_mask.nii \
                  --float \
                  -r ${OUTDIR}/ref_func_bf.nii \
                  -o ${OUTDIR}/wm_mask_nat.nii \
                  -t [${OUTDIR}/func2anat.mat,1] \
                  -t [${OUTDIR}/func2anat0InverseWarp.nii.gz,0] \
                  -n Linear


            ${ABIN}/antsApplyTransforms -v 1 \
                  -i ${OUTDIR}/gm_mask.nii \
                  --float \
                  -r ${OUTDIR}/ref_func_bf.nii \
                  -o ${OUTDIR}/gm_mask_nat.nii \
                  -t [${OUTDIR}/func2anat.mat,1] \
                  -t [${OUTDIR}/func2anat0InverseWarp.nii.gz,0] \
                  -n Linear
      else

            ${ABIN}/antsApplyTransforms -v 1 \
                  -i ${OUTDIR}/csf_mask.nii \
                  --float \
                  -r ${OUTDIR}/ref_func_bf.nii \
                  -o ${OUTDIR}/csf_mask_nat.nii \
                  -t [${OUTDIR}/func2anat.mat,1] \
                  -n Linear

            ${ABIN}/antsApplyTransforms -v 1 \
                  -i ${OUTDIR}/wm_mask.nii \
                  --float \
                  -r ${OUTDIR}/ref_func_bf.nii \
                  -o ${OUTDIR}/wm_mask_nat.nii \
                  -t [${OUTDIR}/func2anat.mat,1] \
                  -n Linear

            ${ABIN}/antsApplyTransforms -v 1 \
                  -i ${OUTDIR}/gm_mask.nii \
                  --float \
                  -r ${OUTDIR}/ref_func_bf.nii \
                  -o ${OUTDIR}/gm_mask_nat.nii \
                  -t [${OUTDIR}/func2anat.mat,1] \
                  -n Linear
      fi

      # ERODE tissue masks
      fslmaths ${OUTDIR}/csf_mask_nat -thr 0.66 -bin  ${OUTDIR}/csf_mask_nat_ero
      gunzip -f ${OUTDIR}/csf_mask_nat_ero.nii.gz

      fslmaths ${OUTDIR}/wm_mask_nat  -bin -kernel 3D -ero ${OUTDIR}/wm_mask_nat_ero
      gunzip -f ${OUTDIR}/wm_mask_nat_ero.nii.gz

      3dcalc -a ${OUTDIR}/wm_mask_nat_ero.nii -b ${OUTDIR}/csf_mask_nat_ero.nii -expr 'a+b' -prefix ${OUTDIR}/nongm_mask_ero.nii

      fslmeants -i ${OUTDIR}/proc_data_native.nii -m ${OUTDIR}/csf_mask_nat_ero.nii -o ${OUTDIR}/csf_sig.txt
      fslmeants -i ${OUTDIR}/proc_data_native.nii -m ${OUTDIR}/wm_mask_nat_ero.nii -o ${OUTDIR}/wm_sig.txt
      fslmeants -i ${OUTDIR}/proc_data_native.nii -m ${OUTDIR}/nat_mask.nii -o ${OUTDIR}/global_sig.txt


      python run_acompcor.py -d ${OUTDIR} -i proc_data_native.nii -n nongm_mask_ero.nii -b nat_mask.nii -t ${TR}

      # Remove headers from compcor results
      sed '1d' ${OUTDIR}/acompcor.txt > ${OUTDIR}/tmp_acompcor.txt
      mv ${OUTDIR}/tmp_acompcor.txt ${OUTDIR}/acompcor.txt

      sed '1d' ${OUTDIR}/tcompcor.txt > ${OUTDIR}/tmp_tcompcor.txt
      mv ${OUTDIR}/tmp_tcompcor.txt ${OUTDIR}/tcompcor.txt

      python rp_nuis_calc.py -d ${OUTDIR} -r motion_estimate.par -t ${FD_THR}

      } &> ${OUTDIR}/03_SignalExtraction.log

      END=$(date -u +%s.%N)
      DIFF=`echo "( $END - $START )" | bc`
      printf "DONE [%.1f s]\n" $DIFF

fi

# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------



# From Smith et al. NeuroImage 2013, HCP

# Preliminary analyses of measures of motion-related artefacts indicate that the ICA-FIX process greatly reduces but does not, in some datasets, totally eliminate motion artefacts
# that are frame-specific and non-spatially-specific.
# In coming months we will investigate further whether there is value in additional cleanup stages,
# most likely to be added into the pipeline after the ICA+FIX denoising. One approach that is simple but effective is “motion scrubbing”, #
# in which one identifies the timepoints that are “irreversibly” damaged by motion, and simply excises those from the timeseries analysis (Power et al., 2011).
#  We will evaluate this and other approaches, and where appropriate, make improvements to the temporal pre-processing pipeline.
#
# Besides this, Dicer [Aquino et al. 2019] supposedly helps with those more global artifacts

if [ "$DO_ICA" -eq "1" ]; then
      printf "[$SUB] Performing MELODIC ICA ... "
      START=$(date -u +%s.%N)

      {

      DOIC=1
      if [ "$DOIC" -eq "1" ]; then
            # Reorder columns in the motion estimates generated by AFNI to FSL's ordering
            rm -f ${OUTDIR}/motion_estimate_fsl.par
            touch ${OUTDIR}/motion_estimate_fsl.par
            while IFS=" " read -r roll pitch yaw ds dl dy
            do
                  rollRad=`echo ${roll} "*3.14159/180" | bc -l | awk '{printf "%f", $0}'`
                  pitchRad=`echo ${pitch} "*3.14159/180" | bc -l | awk '{printf "%f", $0}'`
                  yawRad=`echo ${yaw} "*3.14159/180" | bc -l | awk '{printf "%f", $0}'`

                  echo $rollRad ' ' $pitchRad ' ' $yawRad ' ' $ds ' ' $dl ' ' $dy >> ${OUTDIR}/motion_estimate_fsl.par

            done < ${OUTDIR}/motion_estimate.par

            rm ${OUTDIR}/tmp_melodic.nii
            # Create a temporary detrended 3d+time series to run melodic
            # For ICA-FIX, smoothing is not recommended [only during visualization]
            3dTproject -prefix ${OUTDIR}/tmp_melodic.nii -polort 1 -stopband 0 0.009 -mask ${OUTDIR}/nat_mask.nii  -TR ${TR} -input ${OUTDIR}/proc_data_native.nii

            if test -f "${OUTDIR}/func2anat0Warp.nii.gz"; then
                  ${ABIN}/antsApplyTransforms -v 1 \
                        -i ${OUTDIR}/anat_ref.nii \
                        --float \
                        -r ${OUTDIR}/ref_func_bf.nii \
                        -o ${OUTDIR}/t1_nat_bg.nii \
                        -t [${OUTDIR}/func2anat.mat,1] \
                        -t [${OUTDIR}/func2anat0InverseWarp.nii.gz,0] \
                        -n Linear
            else
                  ${ABIN}/antsApplyTransforms -v 1 \
                        -i ${OUTDIR}/anat_ref.nii \
                        --float \
                        -r ${OUTDIR}/ref_func_bf.nii \
                        -o ${OUTDIR}/t1_nat_bg.nii \
                        -t [${OUTDIR}/func2anat.mat,1] \
                        -n Linear
            fi

            rm -f -R ${OUTDIR}/melodic.ic
            rm -f -R ${OUTDIR}/FIX


            3dpc -eigonly -vmean -vnorm -prefix ${OUTDIR}/pc_var ${OUTDIR}/tmp_melodic.nii
            IDX=0
            NIC=1
            PREV=-1
            while IFS=" " read -r value1 value2 value3 value4
            do
                  if (( IDX > 0 )); then
                        if (( $(echo "${value4} > 0.8" | bc -l) )); then
                          PCT=$(echo "${value4} * 100" | bc -l)
                          NIC=${IDX}
                          #echo "$((IDX-1)) components (${PCT}% of variance explained)"
                          break
                        fi

                        if (( $(echo "${PREV} < 0" | bc -l) )); then
                              PREV=${value3}
                        else

                              CHG=$(echo "${PREV} - ${value3}" | bc -l)

                              if (( $(echo "${CHG} < 0.000001" | bc -l) )); then
                                  echo "[$IDX] Change = ${CHG}"
                                  PCT=$(echo "${value4} * 100" | bc -l)
                                  NIC=$IDX
                                  break
                              fi

                              PREV=${value3}
                        fi

                  fi


                  IDX=$((IDX+1))

            done < "${OUTDIR}/pc_var_eig.1D"

            echo "EXTRACTING ${NIC} components (${PCT}% of variance explained) using ${ICA_TYPE} estimation"

            # --dimest=lap is the default, but it seems empirically to overestimate components
            # see also Varoquaux et al. 2010 NeuroImage


            if [ "$ICA_TYPE" == "melodic" ]; then
                  melodic -i ${OUTDIR}/tmp_melodic.nii -o ${OUTDIR}/melodic.ic --tr=${TR} --mmthresh=0.5 --nobet --dim=${NIC} \
                         --mask=${OUTDIR}/nat_mask.nii --Ostats --report  --eps=0.0001 --bgimage=${OUTDIR}/t1_nat_bg.nii
            fi

            if [ "$ICA_TYPE" == "canica" ]; then
                  mkdir ${OUTDIR}/melodic.ic
                  python ${WDIR}/canica_estimation.py -o ${OUTDIR}/melodic.ic -in ${OUTDIR}/tmp_melodic.nii \
                          -mask ${OUTDIR}/nat_mask.nii -nIC ${NIC}

                  mv ${OUTDIR}/melodic.ic/CanICA_Comps_Z.nii ${OUTDIR}/melodic.ic/melodic_IC.nii
                  gzip ${OUTDIR}/melodic.ic/melodic_IC.nii
                  rm ${OUTDIR}/melodic.ic/melodic_IC.nii

                  # Remove NaNs from the IC maps, in case they are there
                  fslmaths ${OUTDIR}/melodic.ic/melodic_IC.nii.gz -nan ${OUTDIR}/melodic.ic/melodic_IC.nii.gz

                  echo "1" > ${OUTDIR}/melodic.ic/grot.txt
                  melodic -i ${OUTDIR}/melodic.ic/melodic_IC.nii.gz --ICs=${OUTDIR}/melodic.ic/melodic_IC.nii \
                      --mix=${OUTDIR}/melodic.ic/grot.txt -o ${OUTDIR}/melodic.ic --Oall --report -v --mmthresh=0.5


                  #TODO JUST FOR TESTING REMOVE REMOVE REMOVE
                  # melodic -i ${OUTDIR}/tmp_melodic.nii -o ${OUTDIR}/melodic_fsl.ic --tr=${TR} --mmthresh=0.5 --nobet --dim=${NIC} \
                  #      --mask=${OUTDIR}/nat_mask.nii --Ostats --report  --eps=0.0001

                  #3dmerge -doall -1blur_fwhm 5 -prefix ${OUTDIR}/melodic_fsl.ic/melodic_IC_smooth.nii ${OUTDIR}/melodic_fsl.ic/melodic_IC.nii.gz
            fi


            rm ${OUTDIR}/tmp_melodic.nii

            # Creates a smoothed copy of the melodic components. For visualization purposes only
            3dmerge -doall -1blur_fwhm 5 -prefix ${OUTDIR}/melodic.ic/melodic_IC_smooth.nii ${OUTDIR}/melodic.ic/melodic_IC.nii.gz


            # Creates a folder to be used with FIX
            # TODO avoid creating this if ICA-FIX procedure is not going to be used
            mkdir ${OUTDIR}/FIX
            mkdir ${OUTDIR}/FIX/mc
            mkdir ${OUTDIR}/FIX/reg
            mkdir ${OUTDIR}/FIX/filtered_func_data.ica



            cp ${OUTDIR}/FIX/filtered_func_data.ica/melodic_IC.nii.gz ${OUTDIR}/FIX/filtered_func_data.ica/melodic_IC_correct.nii.gz
            fslorient -copysform2qform ${OUTDIR}/FIX/filtered_func_data.ica/melodic_IC_correct.nii.gz
            fslorient -forceneurological ${OUTDIR}/FIX/filtered_func_data.ica/melodic_IC_correct.nii.gz

            fslorient -getorient ${OUTDIR}/FIX/filtered_func_data.ica/melodic_IC.nii.gz
            fslorient -getorient ${OUTDIR}/FIX/filtered_func_data.ica/CanICA_Comps.nii
            fslorient -getorient ${OUTDIR}/FIX/filtered_func_data.ica/melodic_IC_correct.nii.gz

            mv ${OUTDIR}/FIX/filtered_func_data.ica/melodic_IC_correct.nii.gz ${OUTDIR}/FIX/filtered_func_data.ica/melodic_IC.nii.gz


            # This step ensures that func2anat registation mat is properly converted to FSL space
            ${C3DPATH}/c3d_affine_tool -ref ${OUTDIR}/anat_ref.nii -src ${OUTDIR}/ref_func_bfm.nii \
                                -itk ${OUTDIR}/func2anat.mat -ras2fsl -o ${OUTDIR}/afunc2anat_fsl.mat
            convert_xfm -omat ${OUTDIR}/iafunc2anat.mat -inverse ${OUTDIR}/afunc2anat_fsl.mat

            cp ${OUTDIR}/iafunc2anat.mat ${OUTDIR}/FIX/reg/highres2example_func.mat

            cp ${OUTDIR}/proc_data_native.nii ${OUTDIR}/FIX/filtered_func_data.nii
            gzip ${OUTDIR}/FIX/filtered_func_data.nii
            cp ${OUTDIR}/motion_estimate_fsl.par ${OUTDIR}/FIX/mc/prefiltered_func_data_mcf.par
            cp -R ${OUTDIR}/melodic.ic/* ${OUTDIR}/FIX/filtered_func_data.ica
            cp ${OUTDIR}/nat_mask.nii ${OUTDIR}/FIX/mask.nii
            gzip ${OUTDIR}/FIX/mask.nii
            fslmaths ${OUTDIR}/proc_data_native.nii -Tmean ${OUTDIR}/FIX/mean_func
            cp ${OUTDIR}/FIX/mean_func.nii.gz ${OUTDIR}/FIX/reg/example_func.nii.gz
            cp ${OUTDIR}/t1_frn.nii ${OUTDIR}/FIX/reg/highres.nii
            gzip ${OUTDIR}/FIX/reg/highres.nii

            #${C3DPATH}/c3d_affine_tool -ref ${OUTDIR}/ref_func_bf.nii -src ${OUTDIR}/t1_frn.nii ${OUTDIR}/func2anat.mat -inv -oitk ${OUTDIR}/FIX/highres2example_func.mat
            #${C3DPATH}/c3d_affine_tool ${OUTDIR}/afunc2anat.mat -inv -o ${OUTDIR}/iafunc2anat.mat
            #cp ${OUTDIR}/iafunc2anat.mat ${OUTDIR}/FIX/reg/highres2example_func.mat

            CURDIR=`pwd`
            cd ${FIXBIN}
            source ${FIXBIN}/fix -f ${OUTDIR}/FIX
            cd ${CURDIR}

            rm -r -f  ${OUTDIR}/melodic.ic
      fi


      } &> ${OUTDIR}/04_ICA.log

      END=$(date -u +%s.%N)
      DIFF=`echo "( $END - $START )" | bc`
      printf "DONE [%.1f s]\n" $DIFF

fi



# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------

# Use a trained classifier to reduce noise on functional images
# Obviously, this needs a previously trained classifier (and processed func data)
if [ "$APPLY_FIX" -eq "1" ]; then
        printf "[$SUB] Applying FIX Classifier ... "
        START=$(date -u +%s.%N)
        {
        rm -f ${OUTDIR}/proc_data_native_fix.nii
        rm -f ${OUTDIR}/proc_data_native_fixp.nii

        CURDIR=`pwd`
        cd ${FIXBIN}
        source ${FIXBIN}/fix -c ${OUTDIR}/FIX ${FIX_CLASSIFIER} ${FIX_THR}

        # Non-aggressive clean-up
        source ${FIXBIN}/fix -a ${OUTDIR}/FIX/fix4melview_${FIX_CL_LABEL}_thr${FIX_THR}.txt
        cd ${CURDIR}

        mv ${OUTDIR}/FIX/filtered_func_data_clean.nii.gz ${OUTDIR}/proc_data_native_fix.nii.gz
        gunzip -f ${OUTDIR}/proc_data_native_fix.nii.gz
        python ${QCDIR}/save_image_diff.py -o ${OUTDIR} -i ${OUTDIR} -a proc_data_native.nii -b  proc_data_native_fix.nii -type std -msg1 'Before FIX' -msg2 'After FIX'

        # Physiological noise estimation using PHYCAA+ [TODO REF]
        if [ "${DO_PEST}" -eq "2" ]; then
              matlab  "-nodesktop -nosplash " <<<"apply_phycaa_rsn(${TR}, '${OUTDIR}', 'proc_data_native_fix', '${OUTDIR}/nat_mask.nii', '${OUTDIR}/csf_mask_nat.nii' ); exit;"
        fi
      }  &> ${OUTDIR}/07_FIX_apply.log
      END=$(date -u +%s.%N)
      DIFF=`echo "( $END - $START )" | bc`
      printf "DONE [%.1f s]\n" $DIFF
fi


# APPLY_DICER = 1 --> Apply DiCER to native func images
if [ "$APPLY_DICER" -eq "1" ]; then
        printf "[$SUB] DiCER ... "
        START=$(date -u +%s.%N)
        {
        # DICER
        # See how to best do this
        # In the VM, 2mm^3 voxel seems to be too much in terms of memory requirements...
        # Write a regressor for selected ICs, might not be needed
        CLASS_FILE="${OUTDIR}/FIX/fix4melview_${FIX_CL_LABEL}_thr${FIX_THR}.txt"
        COMP_LIST=`tail -1 ${CLASS_FILE}`
        COMP_LIST=${COMP_LIST#"["}

        COMP_LIST=${COMP_LIST%"]"}

        rm ${OUTDIR}/ic_tcs.txt

        IDX=0
        IFS=',' read -a COMPS <<< "$COMP_LIST"
        for COMP in "${COMPS[@]}"; do
            COMP=${COMP//[[:blank:]]/}

            rm ${OUTDIR}/tmp.txt
            touch ${OUTDIR}/tmp.txt

            FSTRING=`tail "${OUTDIR}/FIX/filtered_func_data.ica/report/t${COMP}.txt"`
            while IFS=" " read -r line
            do
                IFS=' ' read -a ELS <<< "$line"
                for EL in "${ELS[@]}"; do
                    echo $EL >> ${OUTDIR}/tmp.txt
                done
                #/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/RS001/FIX/filtered_func_data.ica/report/t2.txt
            done <<<${FSTRING}

            if [ "${IDX}" -eq "0" ]; then
                  IDX=$((IDX+1))
                  cp ${OUTDIR}/tmp.txt ${OUTDIR}/ic_tcs.txt
            else
                  IDX=$((IDX+1))
                  paste  ${OUTDIR}/ic_tcs.txt ${OUTDIR}/tmp.txt > ${OUTDIR}/ic_tcs_tmp.txt
                  mv ${OUTDIR}/ic_tcs_tmp.txt ${OUTDIR}/ic_tcs.txt
            fi
          done


        python rp_nuis_calc.py -d ${OUTDIR} -r motion_estimate.par -t ${FD_THR}
        ${ABIN}/antsApplyTransforms -v 0 \
                    -i ${OUTDIR}/t1_frn.nii \
                    --float \
                    -r ${OUTDIR}/ref_func_bfm.nii \
                    -o ${OUTDIR}/anat2group_lowres.nii \
                    -t [${OUTDIR}/anat2group1Warp.nii.gz,0] \
                    -t [${OUTDIR}/anat2group0GenericAffine.mat,0] \
                    -n BSpline

        if [ "${DO_ANAT_FMAP}" == "1" ]; then

            ${ABIN}/antsApplyTransforms -v 0 \
                      -i ${OUTDIR}/t1_frn.nii \
                      --float \
                      -r ${OUTDIR}/ref_func_bf.nii \
                      -o ${OUTDIR}/anat2group_lowres.nii \
                      -t [${OUTDIR}/func2anat.mat,1] \
                      -t [${OUTDIR}/func2anat0InverseWarp.nii.gz,0] \
                      -n BSpline

        else
            ${ABIN}/antsApplyTransforms -v 0 \
                      -i ${OUTDIR}/t1_frn.nii \
                      --float \
                      -r ${OUTDIR}/ref_func_bf.nii \
                      -o ${OUTDIR}/anat2group_lowres.nii \
                      -t [${OUTDIR}/func2anat.mat,1] \
                      -n BSpline
        fi




        rm -f -r ${OUTDIR}/dicer
        rm -f -r ${OUTDIR}/dicer_fix
        CURDIR=`pwd`

        # TODO Set up conf for DiCER PATH
        cd ${DICERPATH}


        mkdir ${OUTDIR}/dicer_fix
        #mkdir ${OUTDIR}/dicer

        cp ${OUTDIR}/motion_regressors_12_24.txt ${OUTDIR}/dicer/regressors.txt
        cp ${OUTDIR}/motion_regressors_12_24.txt ${OUTDIR}/dicerregressors.txt

        # If background is not 0, fast has problem getting the correct tissue distrub
        # TODO set this value based on quantile
        fslmaths ${OUTDIR}/anat2group_lowres.nii -thr 10 ${OUTDIR}/anat2group_lowres_thr.nii

        #3dcalc -a ${OUTDIR}/proc_data_native_fix.nii -b ${OUTDIR}/gm_mask_nat.nii -expr "a*b" -prefix ${OUTDIR}/tmp.nii
        cp ${OUTDIR}/proc_data_native_fix.nii ${OUTDIR}/tmp.nii

        source DiCER_lightweight.sh -i ${OUTDIR}/tmp -a ${OUTDIR}/anat2group_lowres_thr.nii.gz -w ${OUTDIR}/dicer_fix -s SUBJECT_1_FIX  -p ${DOWN} -d
        #source DiCER_lightweight.sh -i ${OUTDIR}/proc_data_native -a ${OUTDIR}/anat2group_lowres_thr.nii.gz -w ${OUTDIR}/dicer -s SUBJECT_1  -p 2 -d -c regressors.txt
        cd ${CURDIR}

        mv ${OUTDIR}/dicer_fix/tmp_detrended_hpf_dbscan.nii.gz  ${OUTDIR}/proc_data_native_fix_d.nii.gz
        #mv ${OUTDIR}/dicer_fix/proc_data_native_fix_detrended_hpf_dbscan.nii.gz  ${OUTDIR}/proc_data_native_fix_d.nii.gz
        #mv ${OUTDIR}/dicer/proc_data_native_detrended_hpf_dbscan.nii.gz  ${OUTDIR}/proc_data_native_fix_d.nii.gz

        gunzip -f ${OUTDIR}/proc_data_native_fix_d.nii.gz
        #gunzip -f ${OUTDIR}/proc_data_native_d.nii.gz


        rm -f -r ${OUTDIR}/dicer
        rm -f -r ${OUTDIR}/dicer_fix
        rm -f ${OUTDIR}/tmp.nii
        python ${QCDIR}/save_image_diff.py -o ${OUTDIR} -i ${OUTDIR} -a proc_data_native_fix.nii -b  proc_data_native_fix_d.nii -type std -msg1 'Before DiCER' -msg2 'After DiCER'
        #python ${QCDIR}/save_image_diff.py -o ${OUTDIR} -i ${OUTDIR} -a proc_data_native.nii -b  proc_data_native_d.nii -type std -msg1 'Before DiCER' -msg2 'After DiCER'
      }  &> ${OUTDIR}/08_DiCER.log
      END=$(date -u +%s.%N)
      DIFF=`echo "( $END - $START )" | bc`
      printf "DONE [%.1f s]\n" $DIFF
fi






# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------


# From Parkes et al. 2018
#Taken together, the above results allow us to draw six key conclusions.
# First, Head Motion Parameters + Phys models without GSR are ineffective at
#      mitigating motion-related artefact from rs-fMRI,
#      regardless of the level of motion in the sample,
#      exclusionary criteria applied, or the use of expansion terms.
# The application of GSR dramatically improves the performance of these pipelines.
# Second, aCompCor pipelines may only be viable in low-motion datasets and perform poorly in high motion data (see Fig. 1).
# Third, ICA-AROMA and censoring pipelines are superior to the other denoising strategies,
#      yielding the lowest QC-FC correlations (see Fig. 1, Fig. 4), lowest QC-FC distance-dependence (see Fig. 2, Fig. 4),
#      and minimal functional connectivity differences between high- and low-motion healthy controls (see Fig. 6).
# Fourth, a major part of the benefit to motion-related artefact control in censoring pipelines comes
#      from the exclusion of participants with <4 min of uncensored data;
#      when this criterion is applied to all pipelines,
#      performance differences are marginal (except for the HMP pipelines without GSR).
# Fifth, aCompCor and censoring pipelines yield high tDOF-loss,
#      marking them as relatively expensive methods for controlling
#      for motion-related artefact in rs-fMRI data.
# Finally, methods that were more effective at denoising were associated
#      with reduced test-retest reliability,
#      suggesting that noise signals in BOLD data are reproducible.

if [ "$DO_QA" -eq "2" ]; then
      # Performs QS in native space
      printf "[$SUB] Generating QC plots ... "
      START=$(date -u +%s.%N)

      {




          fslmaths ${OUTDIR}/proc_data_native.nii -Tmean ${OUTDIR}/func_mean.nii
          gunzip -f ${OUTDIR}/func_mean.nii.gz



          model_qa()
          {
              OUTDIR=${1}
              ANAT_REF=${2}
              TR=${3}
              QCDIR=${4}
              INPREF=${5}
              REG_MODEL=${6}

              mkdir ${OUTDIR}/QA_${REG_MODEL}
              QA_DIR=${OUTDIR}/QA_${REG_MODEL}

              rm -f ${OUTDIR}/${INPREF}_${REG_MODEL}.nii
              rm -f ${OUTDIR}/${INPREF}_${REG_MODEL}.nii.gz
              rm -f ${QA_DIR}/*

              CMD="3dTproject -input ${OUTDIR}/${INPREF}.nii -TR ${TR} -norm "
              CMD="${CMD} -mask ${OUTDIR}/gm_mask_nat.nii  "
              CMD="${CMD} -prefix ${OUTDIR}/${INPREF}_${REG_MODEL}.nii "



              CMD="${CMD} -cenmode NTRP -censor ${OUTDIR}/temporal_mask_fd.txt "

              # Model-specific regressors
              if [ "${REG_MODEL}" == "SFIX_D" ]; then
                  1dBandpass -dt ${TR} 0.009 100 ${OUTDIR}/motion_estimate.par > ${OUTDIR}/motion_estimate_bp.par
                  CMD="${CMD} -ort ${OUTDIR}/motion_estimate_bp.par "
              fi


              if [ "${REG_MODEL}" == "SRP24WM1CSF1" ] || [ "${REG_MODEL}" == "SRP24CC" ] || \
                 [ "${REG_MODEL}" == "SFIX" ]; then
                  CMD="${CMD} -ort ${OUTDIR}/motion_regressors_12.txt "
              fi

              if [ "${REG_MODEL}" == "SFIX_CC" ]; then
                  #fslmeants -i ${OUTDIR}/${INPREF}.nii -m ${OUTDIR}/csf_mask_nat.nii > ${OUTDIR}/csf_${REG_MODEL}.txt
                  #fslmeants -i ${OUTDIR}/${INPREF}.nii -m ${OUTDIR}/wm_mask_nat_ero.nii > ${OUTDIR}/wm_${REG_MODEL}.txt
                  rm -f ${OUTDIR}/dvars.txt
                  #3dTto1D -input ${OUTDIR}/${INPREF}.nii -prefix ${OUTDIR}/dvars.txt -method srms
                  python run_acompcor.py -d ${OUTDIR} -i ${INPREF}.nii -n nongm_mask_ero.nii -b nat_mask.nii -t ${TR} \
                                        -aout acompcor_fix -tout tcompcor_fix -var 0.8

                  # Remove headers from compcor results
                  sed '1d' ${OUTDIR}/acompcor_fix.txt > ${OUTDIR}/tmp_acompcor_fix.txt
                  mv ${OUTDIR}/tmp_acompcor_fix.txt ${OUTDIR}/acompcor_fix.txt

                  sed '1d' ${OUTDIR}/tcompcor_fix.txt > ${OUTDIR}/tmp_tcompcor_fix.txt
                  mv ${OUTDIR}/tmp_tcompcor_fix.txt ${OUTDIR}/tcompcor_fix.txt

                  python fix_acompcor_pc.py -compcor ${OUTDIR}/acompcor_fix.txt \
                      -fix_ics ${OUTDIR}/tcompcor_fix.txt -expVar 0.25 -svar 2 \
                      -mot12 ${OUTDIR}/motion_regressors_12.txt \
                      -mot24 ${OUTDIR}/motion_regressors_24.txt \
                      -out ${OUTDIR}/fix_cc_nuis.txt

                  #CMD="${CMD} -ort ${OUTDIR}/ic_tcs.txt "
                  #CMD="${CMD} -ort ${OUTDIR}/motion_estimate.par "
                  CMD="${CMD} -ort ${OUTDIR}/fix_cc_nuis.txt "
                  #CMD="${CMD} -ort ${OUTDIR}/dvars.txt "
              fi

              if [ "${REG_MODEL}" == "SRP9" ] || [ "${REG_MODEL}" == "SRP24WM1CSF1" ] \
                || [ "${REG_MODEL}" == "SFIX" ]; then
                  CMD="${CMD} -ort ${OUTDIR}/csf_sig.txt "
                  CMD="${CMD} -ort ${OUTDIR}/wm_sig.txt "
              fi

              if [ "${REG_MODEL}" == "SRP24WM1CSF1" ] || [ "${REG_MODEL}" == "SRP24CC" ] \
                 [ "${REG_MODEL}" == "SFIX" ]; then
                  CMD="${CMD} -ort ${OUTDIR}/motion_regressors_24.txt "
              fi

              if [ "${REG_MODEL}" == "SRP9" ]; then
                  CMD="${CMD} -ort ${OUTDIR}/global_sig.txt "
                  CMD="${CMD} -ort ${OUTDIR}/motion_estimate.par "
              fi

              if [ "${REG_MODEL}" == "SRP24CC" ]; then
                  CMD="${CMD} -ort ${OUTDIR}/acompcor.txt "
              fi

              if [ "${REG_MODEL}" == "SRP24WM1CSF1" ] || [ "${REG_MODEL}" == "SFIX" ] \
                 [ "${REG_MODEL}" == "SRP24CC" ] || [ "${REG_MODEL}" == "SRP9" ] || \
                 [ "${REG_MODEL}" == "NONE" ] || [ "${REG_MODEL}" == "SFIX_CC" ]; then
                  CMD="${CMD} -stopband 0 0.01 -polort 2 "
              else
                  CMD="${CMD} -stopband 0 0.01 -polort 0 "
              fi

              echo "RUNNING: ${CMD}"
              eval "${CMD}"
              rm -f ${OUTDIR}/csf_${REG_MODEL}.txt
              rm -f ${OUTDIR}/wm_${REG_MODEL}.txt


              #python ${QCDIR}/QC_grey_plot.py -out ${QA_DIR} -in ${OUTDIR}  -mpe motion_estimate.par \
              #                  -fname ${INPREF}_${REG_MODEL}.nii -csf_name csf_mask_nat.nii -wm_name wm_mask_nat.nii -gm_name gm_mask_nat.nii \
              #                  -norm zscore -range 1.0 \
              #                  -prog AFNI -outf 02_Greyplot_${REG_MODEL}  -dpi 100 -tr ${TR}

              #python ${QCDIR}/QC_temporal_stats.py -out ${QA_DIR} -in ${OUTDIR}  -mpe motion_estimate.par \
              #                  -fname ${INPREF}_${REG_MODEL}.nii -bg ${ANAT_REF} -type GlobalCorr \
              #                  -range 95% -plane axial -thr 0.3 -smooth 0 -save_nii 1 \
              #                  -prog AFNI -outf 04_GlobalCorr_z_${REG_MODEL}  -dpi 100

              #python ${QCDIR}/QC_temporal_stats.py -out ${QA_DIR} -in ${OUTDIR}  -mpe1 motion_estimate.par \
              #                  -fname ${INPREF}_${REG_MODEL}.nii -bg ${ANAT_REF} -type MotionCorr \
              #                  -range 95% -plane axial -thr 0.3 -smooth 0 -save_nii 1 \
              #                  -prog AFNI -outf 04_MotionCorr_z_${REG_MODEL}  -dpi 100

              python ${QCDIR}/QC_FC.py -out ${QA_DIR} -in ${OUTDIR}  \
                          -fname ${INPREF}_${REG_MODEL}.nii -outf 05_FC_${REG_MODEL}  \
                          -aal ${OUTDIR}/aal2_nat_atlas.nii \
                          -localGlobal ${OUTDIR}/local_global400_nat_atlas.nii \
                          -vmin -0.6 -vmax 0.6 \
                           -dpi 72

              rm -f ${OUTDIR}/${INPREF}_${REG_MODEL}.nii
              rm -f ${OUTDIR}/${INPREF}_${REG_MODEL}.nii.gz
      }
      export -f model_qa
      #parallel -j1 --line-buffer model_qa ::: ${OUTDIR} ::: ${OUTDIR}/t1_nat_bg.nii ::: ${TR} ::: ${QCDIR}  :::  proc_data_native ::: SRP24WM1CSF1
      parallel -j4 --line-buffer model_qa ::: ${OUTDIR} ::: ${OUTDIR}/t1_nat_bg.nii ::: ${TR} ::: ${QCDIR}  :::  proc_data_native ::: NONE SRP24WM1CSF1 SRP24CC SRP9
      #parallel -j4 --line-buffer model_qa ::: ${OUTDIR} ::: ${OUTDIR}/t1_nat_bg.nii ::: ${TR} ::: ${QCDIR}  :::  proc_data_native_d ::: SRP24WM1CSF1_D SRP24_D SFIX_D
      #parallel -j2 --line-buffer model_qa ::: ${OUTDIR} ::: ${OUTDIR}/t1_nat_bg.nii ::: ${TR} ::: ${QCDIR}  :::  proc_data_native_fix ::: SFIX SFIX_CC
      parallel -j1 --line-buffer model_qa ::: ${OUTDIR} ::: ${OUTDIR}/t1_nat_bg.nii ::: ${TR} ::: ${QCDIR}  :::  proc_data_native_fix ::: SFIX SFIX_CC
      parallel -j1 --line-buffer model_qa ::: ${OUTDIR} ::: ${OUTDIR}/t1_nat_bg.nii ::: ${TR} ::: ${QCDIR}  :::  proc_data_native_fix_d ::: SFIX_D


} &> ${OUTDIR}/06_qa.log

      END=$(date -u +%s.%N)
      DIFF=`echo "( $END - $START )" | bc`
      printf "DONE [%.1f s]\n" $DIFF
fi





# Apply normalisation WARP
if [ "$DO_NORM" -eq "1" ]; then

      printf "[$SUB] Warping functional images to MNI space ... "
      START=$(date -u +%s.%N)

      {
            # Run selected nuisance regression
            INPREF=proc_data_native
            if [ "${REG_MODEL}" == "SFIX_CC" ] || [ "${REG_MODEL}" == "SFIX" ]; then
              INPREF=proc_data_native_fix
            fi
            if [ "${REG_MODEL}" == "SFIX_D" ]; then
              INPREF=proc_data_native_fix_d
            fi

            CMD="3dTproject -input ${OUTDIR}/${INPREF}.nii -TR ${TR} "
            CMD="${CMD} -mask ${OUTDIR}/nat_mask.nii "
            CMD="${CMD} -prefix ${OUTDIR}/tmp.nii "



            CMD="${CMD} -cenmode NTRP -censor ${OUTDIR}/temporal_mask_fd.txt "

            # Model-specific regressors
            if [ "${REG_MODEL}" == "SRP24WM1CSF1" ] || [ "${REG_MODEL}" == "SRP24CC" ] || \
               [ "${REG_MODEL}" == "SFIX" ]; then

                CMD="${CMD} -ort ${OUTDIR}/motion_regressors_12.txt "
            fi

            if [ "${REG_MODEL}" == "SFIX_D" ]; then
                1dBandpass -dt ${TR} 0.009 100 ${OUTDIR}/motion_estimate.par > ${OUTDIR}/motion_estimate_bp.par
                CMD="${CMD} -ort ${OUTDIR}/motion_estimate_bp.par "
            fi

            if [ "${REG_MODEL}" == "SFIX_CC" ]; then
                #fslmeants -i ${OUTDIR}/${INPREF}.nii -m ${OUTDIR}/csf_mask_nat.nii > ${OUTDIR}/csf_${REG_MODEL}.txt
                #fslmeants -i ${OUTDIR}/${INPREF}.nii -m ${OUTDIR}/wm_mask_nat_ero.nii > ${OUTDIR}/wm_${REG_MODEL}.txt
                rm -f ${OUTDIR}/dvars.txt
                #3dTto1D -input ${OUTDIR}/${INPREF}.nii -prefix ${OUTDIR}/dvars.txt -method srms
                python run_acompcor.py -d ${OUTDIR} -i ${INPREF}.nii -n nongm_mask_ero.nii -b nat_mask.nii -t ${TR} \
                                      -aout acompcor_fix -tout tcompcor_fix -var 0.5

                # Remove headers from compcor results
                sed '1d' ${OUTDIR}/acompcor_fix.txt > ${OUTDIR}/tmp_acompcor_fix.txt
                mv ${OUTDIR}/tmp_acompcor_fix.txt ${OUTDIR}/acompcor_fix.txt

                sed '1d' ${OUTDIR}/tcompcor_fix.txt > ${OUTDIR}/tmp_tcompcor_fix.txt
                mv ${OUTDIR}/tmp_tcompcor_fix.txt ${OUTDIR}/tcompcor_fix.txt

                python fix_acompcor_pc.py -compcor ${OUTDIR}/acompcor_fix.txt \
                    -fix_ics ${OUTDIR}/tcompcor_fix.txt -expVar 0.33 -svar 2 \
                    -mot12 ${OUTDIR}/motion_regressors_12.txt \
                    -mot24 ${OUTDIR}/motion_regressors_24.txt \
                    -out ${OUTDIR}/fix_cc_nuis.txt

                #CMD="${CMD} -ort ${OUTDIR}/ic_tcs.txt "
                #CMD="${CMD} -ort ${OUTDIR}/motion_estimate.par "
                CMD="${CMD} -ort ${OUTDIR}/acompcor_fix.txt "
                #CMD="${CMD} -ort ${OUTDIR}/dvars.txt "
            fi

            if [ "${REG_MODEL}" == "SRP9" ] || [ "${REG_MODEL}" == "SRP24WM1CSF1" ]; then
                CMD="${CMD} -ort ${OUTDIR}/csf_sig.txt "
                CMD="${CMD} -ort ${OUTDIR}/wm_sig.txt "
            fi

            if [ "${REG_MODEL}" == "SRP24WM1CSF1" ] || [ "${REG_MODEL}" == "SRP24CC" ] \
                || [ "${REG_MODEL}" == "SFIX_D" ]; then
                CMD="${CMD} -ort ${OUTDIR}/motion_regressors_24.txt "
            fi

            if [ "${REG_MODEL}" == "SRP9" ]; then
                CMD="${CMD} -ort ${OUTDIR}/global_sig.txt "
                CMD="${CMD} -ort ${OUTDIR}/motion_estimate.par "
            fi

            if [ "${REG_MODEL}" == "SRP24CC" ]; then
                CMD="${CMD} -ort ${OUTDIR}/acompcor.txt "
            fi

            if [ "${REG_MODEL}" == "SRP24WM1CSF1" ] || [ "${REG_MODEL}" == "SFIX" ] \
               [ "${REG_MODEL}" == "SRP24CC" ] || [ "${REG_MODEL}" == "SRP9" ] || \
               [ "${REG_MODEL}" == "NONE" ] || [ "${REG_MODEL}" == "SFIX_CC" ]; then
                CMD="${CMD} -stopband 0 0.009 -polort 2 "
            else
                CMD="${CMD} -stopband 0 0.009 -polort 0 "
            fi

            echo "RUNNING: ${CMD}"
            #eval "${CMD}"

            rm ${OUTDIR}/${INPREF}_${REG_MODEL}.nii
            rm ${OUTDIR}/proc_data_mni.nii
            #fslmaths ${OUTDIR}/gm_mask_nat.nii -thr 0.25 -bin ${OUTDIR}/gm_mask_nat_bin
            #fslmaths ${OUTDIR}/${INPREF}.nii -Tmean ${OUTDIR}/mean_tmp.nii
            #gunzip -f ${OUTDIR}/mean_tmp.nii
            #${ABIN}/N4BiasFieldCorrection -i ${OUTDIR}/mean_tmp.nii -s 1  -o [ ${OUTDIR}/mean_tmp_b.nii ]

            #3dcalc -a ${OUTDIR}/tmp.nii -b ${OUTDIR}/mean_tmp_b.nii  -expr "a+(b*0.5)" -prefix ${OUTDIR}/${INPREF}_${REG_MODEL}.nii
            #3dcalc -a ${OUTDIR}/proc_data_fix.nii -b ${OUTDIR}/ref_vol.nii -c -expr "a+b" -prefix ${OUTDIR}/${INPREF}_${REG_MODEL}.nii
            #cp ${OUTDIR}/proc_data_native_fix.nii  ${OUTDIR}/${INPREF}_${REG_MODEL}.nii
            #rm -f ${OUTDIR}/tmp.nii
            #rm -f ${OUTDIR}/mean_tmp_b.nii
            #rm -f ${OUTDIR}/mean_tmp.nii

            cp ${OUTDIR}/${INPREF}.nii ${OUTDIR}/${INPREF}_${REG_MODEL}.nii

            norm_func()
            {

                  TN=$1

                  OUTDIR=$2
                  ABIN=$3
                  REF=${4}
                  SMOOTH=${5}
                  CORRECT_ANAT_FMAP_NATIVE=${6}

                  if (( $TN < 10 )); then
                        FI=00$TN
                  elif (( $TN < 100 )); then
                        FI=0$TN
                  else
                        FI=$TN
                  fi

                  if [ "$SMOOTH" -eq "-1" ]; then
                        ${ABIN}/DenoiseImage -x ${OUTDIR}/nat_mask.nii -s 1 -n Rician -i ${OUTDIR}/tmp/proc_data.${FI}.nii  -o [ ${OUTDIR}/tmp/proc_data_s.${FI}.nii ]
                        mv ${OUTDIR}/tmp/proc_data_s.${FI}.nii ${OUTDIR}/tmp/proc_data.${FI}.nii
                  fi

                  if [ "$SMOOTH" -gt "0" ]; then
                        # Smoothing in FSL is defined by sigma [see https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=fsl;54608a1b.1111]
                        #You can use:
                        #fslmaths original.nii -kernel gauss 2.1233226 -fmean smoothed.nii

                        #The gaussian kernel takes its argument as sigma in mm instead of the FWHM, but
                        #you can see how to convert between these two values here:
                        #http://mathworld.wolfram.com/GaussianFunction.html

                        #Basically, divide your FWHM (in mm) by 2.3548 to get sigma.
                        SIGMA=`echo "${SMOOTH}/2.3548" | bc -l`
                        fslmaths ${OUTDIR}/tmp/proc_data.${FI}.nii -kernel gauss ${SIGMA} -fmean ${OUTDIR}/tmp/proc_data_s.${FI}.nii

                        # NOTE: It is better to set FSL to save nii files by default (in terms of speed), though this might create problems for FEAT/FLAME (and certainly for FIX)
                        # I have not tested FSL's pipeline after changing /usr/local/fsl/etc/fslconf/fsl.sh and /usr/local/fsl/etc/fslconf/fsl.csh [OUTPUT type]
                        # Likely, this can be done at the start of this script, which is a TODO
                        gunzip -f ${OUTDIR}/tmp/proc_data_s.${FI}.nii.gz
                        mv ${OUTDIR}/tmp/proc_data_s.${FI}.nii ${OUTDIR}/tmp/proc_data.${FI}.nii
                  fi



                  if [ "${CORRECT_ANAT_FMAP_NATIVE}" -eq "0" ] && test -f "${OUTDIR}/func2anat0Warp.nii.gz"; then
                        ${ABIN}/antsApplyTransforms -v 0 \
                              -i ${OUTDIR}/tmp/proc_data.$FI.nii \
                              --float \
                              -r ${REF} \
                              -o ${OUTDIR}/tmp/proc_data_MNI_$FI.nii \
                              -t [${OUTDIR}/anat2group1Warp.nii.gz,0] \
                              -t [${OUTDIR}/anat2group0GenericAffine.mat,0] \
                              -t [${OUTDIR}/func2anat0Warp.nii.gz,0] \
                              -t [${OUTDIR}/func2anat.mat,0] \
                              -n BSpline
                  else
                        ${ABIN}/antsApplyTransforms -v 0 \
                              -i ${OUTDIR}/tmp/proc_data.$FI.nii \
                              --float \
                              -r ${REF} \
                              -o ${OUTDIR}/tmp/proc_data_MNI_$FI.nii \
                              -t [${OUTDIR}/anat2group1Warp.nii.gz,0] \
                              -t [${OUTDIR}/anat2group0GenericAffine.mat,0] \
                              -t [${OUTDIR}/func2anat.mat,0] \
                              -n BSpline
                  fi

                  # Remove negative values due to BSpline interpolation (mostly in the background of the images)
                  3dmerge -prefix ${OUTDIR}/tmp/proc_data_MNI_thr_$FI.nii -1noneg ${OUTDIR}/tmp/proc_data_MNI_$FI.nii

                  if [ "$SMOOTH" -eq "-1" ]; then
                        ${ABIN}/DenoiseImage -s 1 -n Rician -i ${OUTDIR}/tmp/proc_data_MNI_thr_$FI.nii  -o [ ${OUTDIR}/tmp/proc_data_MNI_thr_s_$FI.nii ]
                        mv ${OUTDIR}/tmp/proc_data_MNI_thr_s_$FI.nii ${OUTDIR}/tmp/proc_data_MNI_thr_$FI.nii
                  fi

            }
            export -f norm_func

            # TODO add other smoothing options here
            SMOOTH=0
            if [ "$DO_NLM" -eq "2" ]; then
                  print_debug 'Nonlocal Means filtering + Normalization'
                  SMOOTH=-1
            fi



            # In case, for some reason, this was not set before
            NVOLS=`fslnvols ${OUTDIR}/proc_data_native.nii`



            mkdir ${OUTDIR}/tmp
            rm -f ${OUTDIR}/tmp.nii
            rm -f ${OUTDIR}/tmp_1.nii
            rm -f ${OUTDIR}/tmp_2.nii

            3dTsplit4D -prefix ${OUTDIR}/tmp/proc_data.nii -keep_datum ${OUTDIR}/${INPREF}_${REG_MODEL}.nii
            parallel -j8 --line-buffer norm_func ::: $(seq 0 ${NVOLS}) ::: ${OUTDIR} ::: ${ABIN} ::: ${MNI_REF_2mm} ::: $SMOOTH ::: ${CORRECT_ANAT_FMAP_NATIVE}
            rm ${OUTDIR}/proc_data_mni.nii
            3dTcat -prefix ${OUTDIR}/proc_data_mni.nii -tr ${TR} ${OUTDIR}/tmp/proc_data_MNI_thr_*.nii

            fslmaths ${OUTDIR}/proc_data_mni.nii -ing 1000 ${OUTDIR}/tmp_1.nii
            gunzip -f ${OUTDIR}/tmp_1.nii.gz
            #3dDespike -nomask -NEW -localedit -cut 2 2.5 -prefix ${OUTDIR}/tmp_2.nii ${OUTDIR}/tmp_1.nii


            #3dmerge -prefix ${OUTDIR}/tmp.nii -1blur_fwhm 8 -1noneg -doall ${OUTDIR}/tmp_2.nii
            #SIGMA=`echo "8/2.3548" | bc -l`
            #fslmaths ${OUTDIR}/tmp_2.nii -kernel gauss ${SIGMA} -fmean ${OUTDIR}/tmp.nii
            #gunzip -f ${OUTDIR}/tmp.nii.gz
            mv ${OUTDIR}/tmp_1.nii ${OUTDIR}/proc_data_mni.nii

            rm -r ${OUTDIR}/tmp
            rm -f ${OUTDIR}/tmp.nii
            rm -f ${OUTDIR}/tmp_1.nii
            rm -f ${OUTDIR}/tmp_2.nii

            ## SMA
            #3dmaskave -q -nball 0 -10 60 6 ${OUTDIR}/proc_data_mni.nii > ${OUTDIR}/sma.txt

            ## L Putamen
            #3dmaskave -q -nball -20 4 10 6 ${OUTDIR}/proc_data_mni.nii > ${OUTDIR}/l_putamen.txt

            ## R Putamen
            #3dmaskave -q -nball 20 4 10 6 ${OUTDIR}/proc_data_mni.nii > ${OUTDIR}/r_putamen.txt

            ## R CRB
            #3dmaskave -q -nball 24 -60 -40 6 ${OUTDIR}/proc_data_mni.nii > ${OUTDIR}/r_crb.txt
            ##L CRB
            #3dmaskave -q -nball -24 -60 -40 6 ${OUTDIR}/proc_data_mni.nii > ${OUTDIR}/l_crb.txt

            fslmaths ${OUTDIR}/proc_data_mni.nii -Tmean -thr 100 -bin ${OUTDIR}/mni_mask

            #rm -f  ${OUTDIR}/reho.nii
            #3dReHo -inset ${OUTDIR}/proc_data_mni.nii -prefix ${OUTDIR}/reho.nii -neigh_RAD 2.3 -mask ${OUTDIR}/mni_mask.nii.gz


            #3dTcorr1D -prefix ${OUTDIR}/Corr_R_CRB.nii -mask ${OUTDIR}/mni_mask.nii.gz ${OUTDIR}/proc_data_mni.nii ${OUTDIR}/r_crb.txt
            #3dTcorr1D -prefix ${OUTDIR}/Corr_L_CRB.nii -mask ${OUTDIR}/mni_mask.nii.gz ${OUTDIR}/proc_data_mni.nii ${OUTDIR}/l_crb.txt

            #3dTcorr1D -prefix ${OUTDIR}/Corr_R_Putamen.nii -mask ${OUTDIR}/mni_mask.nii.gz ${OUTDIR}/proc_data_mni.nii ${OUTDIR}/r_putamen.txt
            #3dTcorr1D -prefix ${OUTDIR}/Corr_L_Putamen.nii -mask ${OUTDIR}/mni_mask.nii.gz ${OUTDIR}/proc_data_mni.nii ${OUTDIR}/l_putamen.txt

            #3dTcorr1D -prefix ${OUTDIR}/Corr_SMA.nii -mask ${OUTDIR}/mni_mask.nii.gz ${OUTDIR}/proc_data_mni.nii ${OUTDIR}/sma.txt

} &> ${OUTDIR}/05_norm2mni.log

      END=$(date -u +%s.%N)
      DIFF=`echo "( $END - $START )" | bc`
      printf "DONE [%.1f s]\n" $DIFF
fi
exit



#cp -R /mnt/hgfs/CRUNCHOUT/RS${SUB}/* ${OUTDIR}/
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------





END_ALL=$(date -u +%s.%N)
DIFF=`echo "( $END_ALL - $START_ALL )" | bc`
printf "DONE in %.1f seconds\n" $DIFF

exit 1

#=====================================
printf "\t[${NET}] Moving to external drive ... "
START=$(date -u +%s.%N)

{

mkdir ${FINALDIR}
rm ${FINALDIR}/*
} &> /dev/null

mv ${OUTDIR}/* ${FINALDIR}/
rm -r -f ${OUTDIR}


END=$(date -u +%s.%N)
DIFF=`echo "( $END - $START )" | bc`
printf "OK [%.1f s]\n" $DIFF
