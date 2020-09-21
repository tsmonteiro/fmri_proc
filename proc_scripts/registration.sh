#!/bin/bash
# Clear terminal screen
#printf "\033c"

OUTDIR=$1
CONFIG=$2

#TODO Add this to the configuration file
WDIR='/home/luna.kuleuven.be/u0101486/workspace/fmri_proc/'
QCDIR="${WDIR}/QualityControl/"
C3DPATH='/home/luna.kuleuven.be/u0101486/Software/itk/ITK/c3d/bin/'
CMTKPATH='/home/luna.kuleuven.be/u0101486/Software/CMTK/cmtk-3.3.1p1/build/bin/'

# PARSE Configuration file
ABIN=$(awk -F\=  '/^ABIN/{print $2}' "${CONFIG}" | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')

EXTRACT_SURF=$(awk -F\=  '/^EXTRACT_SURF/{print $2}' "${CONFIG}" | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
DO_FUNC_SURF=$(awk -F\=  '/^DO_FUNC_SURF/{print $2}' "${CONFIG}" | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
MNI_REF=$(awk -F \= '/^\<MNI_REF\>/{print $2}' "${CONFIG}" | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
MNI_REF_2mm=$(awk -F\=  '/^\<MNI_REF_2mm\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')

DO_CLEAN=$(awk -F\=  '/^\<DO_CLEAN\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')

DO_ATLAS_NAT=$(awk -F\=  '/^\<DO_ATLAS_NAT\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')


F2A_FUNC=$(awk -F\=  '/^\<F2A_FUNC\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')

DO_ANAT_FMAP=$(awk -F\=  '/^\<DO_ANAT_FMAP\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
CORRECT_ANAT_FMAP_NATIVE=$(awk -F\=  '/^\<CORRECT_ANAT_FMAP_NATIVE\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')

DO_DEOBL=$(awk -F\=  '/^\<DO_DEOBL\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')

BREXT_TYPE=$(awk -F\=  '/^\<BREXT_TYPE\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')

USE_TISSUE_PRIOR=$(awk -F\=  '/^\<USE_TISSUE_PRIOR\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
GM_PRIOR=$(awk -F\=  '/^\<GM_PRIOR\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
WM_PRIOR=$(awk -F\=  '/^\<WM_PRIOR\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
CSF_PRIOR=$(awk -F\=  '/^\<CSF_PRIOR\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')


print_time()
{
      HOUR=`date +%T`
      DAY=`date +%F`
      TYPE=$1

      echo "[$TYPE] (${DAY} - ${HOUR} )"
}

if [ "$DO_DEOBL" -eq "1" ]; then
      print_time 'INFO'
      echo "Deobliquing T1 volume (3dresample + 3drefit)."
      {
          3dresample -prefix ${OUTDIR}/t1_f.nii -orient lpi -input ${OUTDIR}/t1.nii
          3drefit -deoblique ${OUTDIR}/t1_f.nii
      } &> /dev/null
else
      cp ${OUTDIR}/t1.nii ${OUTDIR}/t1_f.nii
fi



# TODO Add the possibility to import this from somewhere
# This would be useful when there are multiple functional blocks (no need process anatomical and register to MNI again in those cases)
if [ "${EXTRACT_SURF}" -eq "1" ]; then
    print_time 'INFO'
    echo "Extracting cortical surface with CAT12. [${OUTDIR}/logs/surface_extract.log]"
    {
        matlab "-nodesktop -nosplash " <<<"cd '${WDIR}'; anat_proc_surf('${OUTDIR}', 't1_f.nii'); exit;"
    } &> ${OUTDIR}/logs/surface_extract.log
fi

print_time 'INFO'
echo "Running bias field correction."
{
    ${ABIN}/N4BiasFieldCorrection -i ${OUTDIR}/t1_f.nii -s 4  -o [ ${OUTDIR}/t1_fb.nii,${OUTDIR}/biasfield.nii ]
} &> /dev/null


if [ "$BREXT_TYPE" -eq "1" ]; then
      print_time 'INFO'
      echo "Extracting anatomical brain mask using deepbrain. [${OUTDIR}/logs/brain_extract.log]"
      echo "Please check the quality of ${OUTDIR}/AnatMask.nii."
      {
            # [NOTE]
            # deepbrain (used here) requires an older version of tensorflow (v1x)
            # the way around this is to edit lines 18, 19, and 20 (all with tf) in file
            # /home/luna.kuleuven.be/u0101486/Software/anaconda3/lib/python3.7/site-packages/deepbrain/extractor.py
            # "According to TF 1:1 Symbols Map, in TF 2.0 you should use tf.compat.v1.Session() instead of tf.Session()"
            # Source: https://stackoverflow.com/questions/55142951/tensorflow-2-0-attributeerror-module-tensorflow-has-no-attribute-session
            python ${WDIR}/extract_brain.py -i "${OUTDIR}/t1_fb.nii" -o "${OUTDIR}/AnatMask.nii"
      } &> ${OUTDIR}/logs/brain_extract.log

      if [ ! -f "${OUTDIR}/AnatMask.nii" ]; then
          print_time 'ERROR'
          echo "Brain extraction failed."
          echo "Please check ${OUTDIR}/logs/brain_extract.log for more info."
          echo "Execution will abort."
          exit 1
      fi
fi


if [ "$BREXT_TYPE" -eq "2" ]; then
      #TODO Add FSL's bet extraction procedure
      print_time 'WARN/ERROR'
      echo "Brain extraction using bet has not been implemented yet."
      echo "Please set BREXT_TYPE to 1 in configuration file and execute again."
      echo "Execution will abort. "
      exit 1
fi
print_time 'INFO'
echo "Masking anatomical image."
{
    3dcalc -a ${OUTDIR}/t1_fb.nii -b ${OUTDIR}/AnatMask.nii -expr 'a*b' -prefix ${OUTDIR}/t1_brain.nii
} &>/dev/null

print_time 'INFO'
echo "Denoising T1 volume with ANTs DenoiseImage + 3dUnifize to adjust intensities. [${OUTDIR}/logs/brain_denoise.log]"
echo "This will remove Gaussian noise from the image, but preserves edge."
#TODO Make this step optional
{
    ${ABIN}/DenoiseImage -i ${OUTDIR}/t1_brain.nii  -o [ ${OUTDIR}/t1_frn.nii,${OUTDIR}/t1_noise.nii ]
    3dUnifize -prefix ${OUTDIR}/anat_ref.nii ${OUTDIR}/t1_frn.nii
} &> ${OUTDIR}/logs/brain_denoise.log

print_time 'INFO'
echo "Registering T1 volume to Template. [${OUTDIR}/logs/anat2mni.log]"
# Anat2MNI. The results of this transform are also used to obtain priors for segmentation
SRC=${OUTDIR}/anat_ref.nii
{
    ${ABIN}/antsRegistration -d 3 -r [$MNI_REF,$SRC,0] -v 1 \
                -m MI[$MNI_REF,$SRC,1,32] -t translation[0.1] -c [500,5.e-7,20] \
                -s 3vox -f 3 -l 1 -n BSpline \
                -m MI[$MNI_REF,$SRC,1,32,Regular,0.25] -t rigid[0.1] -c [500,5.e-7,20] \
                -s 3vox -f 3 -l 1 -n BSpline \
                -m MI[$MNI_REF,$SRC,1,32,Regular,0.25] -t affine[0.1] -c [500x100x10,5.e-7,10] \
                -s 2x1x0vox -f 3x2x1 -l 1 -n BSpline \
                -m CC[$MNI_REF,$SRC,1,3] -t SyN[0.2,3] -c [20x10,1.e-7,10] \
                -s 1x0vox -f 2x1 -l 1 -n BSpline \
                -o [${OUTDIR}/anat2group,${OUTDIR}/anat2group.nii]
} &> ${OUTDIR}/logs/anat2mni.log


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
fslmaths "${OUTDIR}/proc_data_native.nii" -Tmedian -thr 0 "${OUTDIR}/ref_func_b.nii"
gunzip -f ${OUTDIR}/ref_func_b.nii.gz
${ABIN}/DenoiseImage -i ${OUTDIR}/ref_func_b.nii  -o [ ${OUTDIR}/ref_func_bf.nii]
3dAutomask -apply_prefix ${OUTDIR}/ref_func_bfm.nii ${OUTDIR}/ref_func_bf.nii



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
