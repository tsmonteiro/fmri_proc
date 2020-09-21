#!/bin/bash
# Clear terminal screen
#printf "\033c"

SUB=$1
CONFIG=$2 #/home/fsluser/Documents/rs_proc/conf_files/rep_impact_belgium.conf
#BLOCK=$3
DOWN=$3

#TODO Add this to the configuration file
WDIR=$(awk -F\=  '/^WDIR/{print $2}' "${CONFIG}" | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
QCDIR=$(awk -F\=  '/^QCDIR/{print $2}' "${CONFIG}" | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
DICERPATH=$(awk -F\=  '/^DICERPATH/{print $2}' "${CONFIG}" | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
C3DPATH=$(awk -F\=  '/^C3DPATH/{print $2}' "${CONFIG}" | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
CMTKPATH=$(awk -F\=  '/^CMTKPATH/{print $2}' "${CONFIG}" | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
ABIN=$(awk -F\=  '/^ABIN/{print $2}' "${CONFIG}" | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
FIXBIN=$(awk -F\=  '/^FIXBIN/{print $2}' "${CONFIG}" | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')


EXTRACT_SURF=$(awk -F\=  '/^EXTRACT_SURF/{print $2}' "${CONFIG}" | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
DO_FUNC_SURF=$(awk -F\=  '/^DO_FUNC_SURF/{print $2}' "${CONFIG}" | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')


DARTEL_TEMPLATE_PREF=$(awk -F \= '/^\<DARTEL_TEMPLATE_PREF\>/{print $2}' "${CONFIG}" | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
DARTEL_TEMPLATE_DIR=$(awk -F\=  '/^\<DARTEL_TEMPLATE_DIR\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')

TMPDIR=$(awk -F\=  '/^\<TMPDIR\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
OUTDIR=$(awk -F\=  '/^\<OUTDIR\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
FINALDIR=$(awk -F\=  '/^\<FINALDIR\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')

DO_COPY=$(awk -F\=  '/^\<DO_COPY\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
MOCO=$(awk -F\=  '/^\<MOCO\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
ICA_TYPE=$(awk -F\=  '/^\<ICA_TYPE\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
DO_REG=$(awk -F\=  '/^\<DO_REG\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
DO_ANAT=$(awk -F\=  '/^\<DO_ANAT\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
DO_ICA=$(awk -F\=  '/^\<DO_ICA\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
DO_NORM=$(awk -F\=  '/^\<DO_NORM\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
DO_MNI=$(awk -F\=  '/^\<DO_MNI\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
DO_QA=$(awk -F\=  '/^\<DO_QA\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
RUN_ICA=$(awk -F\=  '/^\<RUN_ICA\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
DO_CLEAN=$(awk -F\=  '/^\<DO_CLEAN\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
APPLY_FIX=$(awk -F\=  '/^\<APPLY_FIX\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
APPLY_DICER=$(awk -F\=  '/^\<APPLY_DICER\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')

DO_ATLAS_NAT=$(awk -F\=  '/^\<DO_ATLAS_NAT\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
TR=$(awk -F\=  '/^\<TR\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
ACQ_TYPE=$(awk -F\=  '/^\<ACQ_TYPE\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
F2A_FUNC=$(awk -F\=  '/^\<F2A_FUNC\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
DO_SLC=$(awk -F\=  '/^\<DO_SLC\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
DO_FMAP=$(awk -F\=  '/^\<DO_FMAP\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
EES=$(awk -F\=  '/^\<EES\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
DO_PEST=$(awk -F\=  '/^\<DO_PEST\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')


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

GROUPDIR=$(awk -F\=  '/^\<GROUPDIR\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')


FINALDIR=`echo ${FINALDIR//@SUB@/${SUB}} | tr -d '"'`
OUTDIR=`echo ${OUTDIR//@SUB@/${SUB}} | tr -d '"'`
TMPDIR=`echo ${TMPDIR} | tr -d '"'`
ABIN=`echo ${ABIN} | tr -d '"'`
ANAT_PROC=`echo ${ANAT_PROC} | tr -d '"'`



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
      echo '                   *****************************************************************************************'
      echo '                   *****************************************************************************************'
      echo ''
      echo "                                                  (${DAY} - ${HOUR} ) $1"
      echo ''
      echo '                   *****************************************************************************************'
      echo '                   *****************************************************************************************'
      echo ''
      echo ''

}


nlm_smooth()
{

      TN=$1
      OUTDIR=$2
      ABIN=$3
      MASK=$4



      if (( $TN < 10 )); then
            FI=00$TN
      elif (( $TN < 100 )); then
            FI=0$TN
      else
            FI=$TN
      fi

      {
          # -r 3x3x2 -p 2x2x1
          ${ABIN}/DenoiseImage -x ${MASK} -r 3 -n Rician -i ${OUTDIR}/tmp/func_data.${FI}.nii  -o [ ${OUTDIR}/tmp/func_data_n.${FI}.nii ]
      } &> /dev/null

}
export -f nlm_smooth

copy_header()
{
    IMBASE=$1
    IMTARGET=$2

    # Copy all relevant orientation and position fields
    QFORM_CODE=`nifti_tool -quiet -disp_hdr -field qform_code -infiles ${IMBASE}`
    SFORM_CODE=`nifti_tool -quiet -disp_hdr -field sform_code -infiles ${IMBASE}`

    QUATERNB=`nifti_tool -quiet -disp_hdr -field quatern_b -infiles ${IMBASE}`
    QUATERNC=`nifti_tool -quiet -disp_hdr -field quatern_c -infiles ${IMBASE}`
    QUATERND=`nifti_tool -quiet -disp_hdr -field quatern_d -infiles ${IMBASE}`

    QOFFX=`nifti_tool -quiet -disp_hdr -field qoffset_x -infiles ${IMBASE}`
    QOFFY=`nifti_tool -quiet -disp_hdr -field qoffset_y -infiles ${IMBASE}`
    QOFFZ=`nifti_tool -quiet -disp_hdr -field qoffset_z -infiles ${IMBASE}`

    SROWX=`nifti_tool -quiet -disp_hdr -field srow_x -infiles ${IMBASE}`
    SROWY=`nifti_tool -quiet -disp_hdr -field srow_y -infiles ${IMBASE}`
    SROWZ=`nifti_tool -quiet -disp_hdr -field srow_z -infiles ${IMBASE}`

    PIXDIM=`nifti_tool -quiet -disp_hdr -field pixdim -infiles ${IMBASE}`

    SLSTART=`nifti_tool -quiet -disp_hdr -field slice_start -infiles ${IMBASE}`
    SLEND=`nifti_tool -quiet -disp_hdr -field slice_end -infiles ${IMBASE}`
    SLCODE=`nifti_tool -quiet -disp_hdr -field slice_code -infiles ${IMBASE}`

    nifti_tool -mod_hdr -overwrite -infiles ${IMTARGET}  \
          -mod_field qform_code ${QFORM_CODE} -mod_field sform_code ${SFORM_CODE} \
          -mod_field quatern_b ${QUATERNB} -mod_field quatern_c ${QUATERNC} \
          -mod_field quatern_d ${QUATERND} -mod_field qoffset_x ${QOFFX} \
          -mod_field qoffset_y ${QOFFY} -mod_field qoffset_z ${QOFFZ} \
          -mod_field srow_x "${SROWX}" -mod_field srow_y "${SROWY}" \
          -mod_field srow_z "${SROWZ}" -mod_field pixdim "${PIXDIM}" \
          -mod_field slice_start ${SLSTART} -mod_field slice_end ${SLEND} \
          -mod_field slice_code ${SLCODE}
}

correct_anat_fmap()
{

    TN=$1
    OUTDIR=$2
    ABIN=$3
    WARP=$4



    if (( $TN < 10 )); then
          FI=00$TN
    elif (( $TN < 100 )); then
          FI=0$TN
    else
          FI=$TN
    fi

    {

      ${ABIN}/antsApplyTransforms \
            -i ${OUTDIR}/tmp/func_data.${FI}.nii \
            --float \
            -r ${OUTDIR}/tmp/func_data.${FI}.nii \
            -o ${OUTDIR}/tmp/func_data_n.${FI}.nii \
            -t [${OUTDIR}/${WARP}0GenericAffine.mat,0] \
            -t [${OUTDIR}/${WARP}1Warp.nii.gz,0] \
            -n Linear
    } &> /dev/null

}
export -f correct_anat_fmap

ANATDIR=${OUTDIR}/anat
IMGDIR=${OUTDIR}/figures
FIXDIR=${OUTDIR}/FIX
LOGDIR=${OUTDIR}/logs


#TODO MOVE this to appropriate sections...


#GROUPDIR='/home/luna.kuleuven.be/u0101486/workspace/data/ConnectEx/Template/'


# =================================================
#            PREPROCESSING START
# =================================================
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

      rm -f -r ${IMGDIR}
      rm -f -r ${ANATDIR}
      rm -f -r ${LOGDIR}

      mkdir ${IMGDIR}
      mkdir ${ANATDIR}
      mkdir ${LOGDIR}

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

            NVOLS=`fslnvols ${OUTDIR}/func_data.nii`
            NZ=`3dinfo -nk ${OUTDIR}/func_data.nii`


            PREF=''
            #if test ! -f "${OUTDIR}/slice_acq.txt"; then
            # File is not present (normally this should happen because data was not imported from DICOM)
            python ./util/write_slice_timing.py -tr ${TR} -nsl ${NZ} -acq ${ACQ_TYPE} -out ${OUTDIR}/slice_acq.txt
            #fi


            # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            # If movement is extreme, this might help. Then again, it is questionable when exactly to do this.
            if [ "$DO_DPK" -eq "1" ]; then
                  log_command_div 'Despike START'
                  3dDespike -nomask -NEW -localedit -cut 2 4 -prefix ${OUTDIR}/d${PREF}func_data.nii ${OUTDIR}/${PREF}func_data.nii
                  python ${QCDIR}/save_image_diff.py -o ${OUTDIR} -i ${OUTDIR} -a ${PREF}func_data.nii -b d${PREF}func_data.nii -msg1 'Before Despike' -msg2 'After Despike'
                  PREF=d${PREF}
                  mv ${OUTDIR}/*.png ${IMGDIR}
                  mv ${OUTDIR}/*.gif ${IMGDIR}
                  log_command_div 'Despike END'
            fi





            # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            #
            # Motion correction
            #
            # There are 2 methods here: 3dvolreg or SPM INRIAlign
            # My tests seem to confirm the idea that slomoco performs better than 3dvolreg
            # This procedure, however, is extremely computationaly intensive (i.e. it takes hours instead of minutes to complete)
            #
            # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            if [ "$MOCO" == "inria" ]; then
                  echo 'Not implemented. Sorry, everything will break.'
            fi


            #mkdir ${OUTDIR}/tmp
            #3dTsplit4D -prefix ${OUTDIR}/tmp/func_data.nii -keep_datum ${OUTDIR}/${PREF}func_data.nii
            #3dAutomask -prefix ${OUTDIR}/nat_mask.nii -dilate 2  ${OUTDIR}/${PREF}func_data.nii
            #parallel -j6 --line-buffer nlm_smooth ::: $(seq 0 ${NVOLS}) ::: ${OUTDIR} ::: ${ABIN} ::: ${OUTDIR}/nat_mask.nii

            #3dTcat -prefix ${OUTDIR}/s${PREF}func_data.nii -tr ${TR} ${OUTDIR}/tmp/func_data_n*.nii
            #rm -f -r ${OUTDIR}/tmp
            #PREF=s${PREF}

            if [ "$MOCO" == "3dvolreg" ]; then
                  log_command_div 'Motion Correction START'

                  3dAutomask -prefix ${OUTDIR}/nat_mask.nii ${OUTDIR}/func_data.nii
                  #${WDIR}/proc_scripts/align_pestica.sh ${OUTDIR} func_data ${TR} 0
                  #python ${QCDIR}/save_image_diff.py -o ${OUTDIR} -i ${OUTDIR} -a ${PREF}func_data.nii -b p${PREF}func_data.nii -msg1 'RAW' -msg2 'PESTICA'

                  #mv ${OUTDIR}/p${PREF}func_data.nii ${OUTDIR}/${PREF}func_data.nii

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
                              -x_thresh 0.01 -zpad 5 -maxite 100 -1Dfile ${OUTDIR}/motion_estimate.par -maxdisp1D ${OUTDIR}/maximum_disp.1d  \
                               ${OUTDIR}/${PREF}func_data.nii"

                  # Motion correction
                  # Performing the realignment to the mean volume increases the processing time quite a lot, without any major benefit
                  # (at least that is what Oakes et al. 2005 say ['Comparison of fMRI motion correction software tools'], though connectivity was not investigated )
                  3dTcat -prefix ${OUTDIR}/ref_vol.nii -TR ${TR} ${OUTDIR}/${PREF}func_data.nii[${REF_VOL}]
                  #3dmerge -prefix ${OUTDIR}/ref_vol_smooth.nii -1blur_fwhm 1 ${OUTDIR}/ref_vol.nii
                  3dAutomask  -apply_prefix ${OUTDIR}/ref_vol_masked.nii ${OUTDIR}/ref_vol.nii
                   #-savedisp ${OUTDIR}/voxel_displ.nii
                  3dvolreg -Fourier -prefix ${OUTDIR}/r${PREF}func_data.nii -base ${OUTDIR}/ref_vol_masked.nii -rot_thresh 0.005 -delta 0.1 \
                        -x_thresh 0.001 -zpad 10 -maxite 60 -1Dfile ${OUTDIR}/motion_estimate.par -maxdisp1D ${OUTDIR}/maximum_disp.1d -savedisp ${OUTDIR}/voxel_displ.nii \
                         ${OUTDIR}/${PREF}func_data.nii

                  #3dREMLfit -input ${OUTDIR}/r${PREF}func_data.nii \
                  #    -dsort ${OUTDIR}/voxel_displ_DX.nii \
                  #    -dsort ${OUTDIR}/voxel_displ_DY.nii \
                  #    -dsort ${OUTDIR}/voxel_displ_DZ.nii \
                  #    -GOFORIT \
                  #    -Oerrts ${OUTDIR}/r2${PREF}func_data.nii #-Rwherr ${OUTDIR}/p${PREF}func_data.nii

                  #fslmaths ${OUTDIR}/${PREF}func_data.nii -Tmean ${OUTDIR}/mean.nii
                  #3dcalc -a ${OUTDIR}/r2${PREF}func_data.nii -b ${OUTDIR}/mean.nii.gz -expr 'a+b' -prefix ${OUTDIR}/rm${PREF}func_data.nii
                  #mv ${OUTDIR}/rm${PREF}func_data.nii ${OUTDIR}/r2${PREF}func_data.nii



                  #python ${QCDIR}/save_image_diff.py -o ${OUTDIR} -i ${OUTDIR} -a r${PREF}func_data.nii -b r2${PREF}func_data.nii -msg1 'Volreg' -msg2 'Volreg DX'
                  #mv ${OUTDIR}/r2${PREF}func_data.nii ${OUTDIR}/r${PREF}func_data.nii

                  1dplot -thick -volreg -png ${OUTDIR}/motion_estimate_0.png -one ${OUTDIR}/motion_estimate.par
                  1dplot -thick -png ${OUTDIR}/maximum_disp_0.png -one ${OUTDIR}/maximum_disp.1d
                  1dplot -thick -png ${OUTDIR}/maximum_disp_delt_0.png -one ${OUTDIR}/maximum_disp.1d_delt

                  mv ${OUTDIR}/*.png ${IMGDIR}
                  mv ${OUTDIR}/*.gif ${IMGDIR}

                  PREF=r${PREF}
                  log_command_div 'Motion Correction END'

            fi

            if [ "${DO_PEST}" == "4" ]; then
                  3dAutomask -prefix ${OUTDIR}/ref_mask.nii -dilate 1 ${OUTDIR}/${PREF}func_data.nii
                  3dmerge -prefix ${OUTDIR}/tmp_smooth.nii -doall -1blur_fwhm 1 ${OUTDIR}/${PREF}func_data.nii


                  matlab  "-nodesktop -nosplash " <<<"apply_phycaa_rsn(${TR}, '${OUTDIR}', 'tmp_smooth', '${OUTDIR}/ref_mask.nii', '', 'physio.txt' ); exit;"

                  rm -f ${OUTDIR}/ptmp_smooth.nii
                  rm -f ${OUTDIR}/tmp_smooth.nii

                  rm -f ${OUTDIR}/p${PREF}func_data.nii
                  #3dREMLfit -input ${OUTDIR}/${PREF}func_data.nii \
                  #    -addbase ${OUTDIR}/physio_regressor.txt \
                  #    -GOFORIT \
                  #    -Oerrts ${OUTDIR}/p${PREF}func_data.nii #-Rwherr ${OUTDIR}/p${PREF}func_data.nii

                  #fslmaths ${OUTDIR}/${PREF}func_data.nii -Tmean ${OUTDIR}/mean.nii
                  #3dcalc -a ${OUTDIR}/${PREF}func_data.nii -b ${OUTDIR}/phyca_data__NN_map.nii -expr 'a*(astep(b, 0.9))' -prefix ${OUTDIR}/p${PREF}func_data.nii
                  #fslmaths ${OUTDIR}/phyca_data__NN_map.nii -binv ${OUTDIR}/physio_mask.nii
                  #fslmeants -i ${OUTDIR}/${PREF}func_data.nii -o ${OUTDIR}/physio.txt -m ${OUTDIR}/physio_mask.nii \
                  #      --eig --order=5 --no_bin
                  #3dAutomask -apply_prefix ${OUTDIR}/mp${PREF}func_data.nii ${OUTDIR}/p${PREF}func_data.nii
                  #mv ${OUTDIR}/mp${PREF}func_data.nii ${OUTDIR}/p${PREF}func_data.nii

                  rm -f ${OUTDIR}/phyca_data__Physio_Zmap.nii
                  rm -f ${OUTDIR}/phyca_data__PHYCAA_step1+2.mat
                  #PREF=p${PREF}

                  mv ${OUTDIR}/*.png ${IMGDIR}
                  mv ${OUTDIR}/*.gif ${IMGDIR}
            fi

            if [ "${DO_PEST}" == "3" ]; then
                  log_command_div 'PHYSIO CORRECTION'
                  # At the moment, only working for one specific data set
                  # 4 --> Number of dummy scans. Pass this to parameters
                  # 607 is scan dependent
                  matlab "-nodesktop -nosplash " <<<"cd '${WDIR}/matlab'; process_physio_data('${OUTDIR}', $NVOLS, 607, ${TR}, ${NZ}, -1); exit;"


                  mv ${OUTDIR}/slice_timing_g.txt ${OUTDIR}/slice_acq.txt
                  #python ./util/write_slice_timing.py -tr ${TR} -nsl ${NZ} -acq 3 -out ${OUTDIR}/slice_acq.txt


                  3dREMLfit -input ${OUTDIR}/${PREF}func_data.nii \
                      -addbase ${OUTDIR}/physio_reg.ortho \
                      -GOFORIT \
                      -Oerrts ${OUTDIR}/p${PREF}func_data.nii #-Rwherr ${OUTDIR}/p${PREF}func_data.nii

                  fslmaths ${OUTDIR}/${PREF}func_data.nii -Tmean ${OUTDIR}/mean.nii
                  3dcalc -a ${OUTDIR}/p${PREF}func_data.nii -b ${OUTDIR}/mean.nii.gz -expr 'a+b' -prefix ${OUTDIR}/pm${PREF}func_data.nii
                  mv ${OUTDIR}/pm${PREF}func_data.nii ${OUTDIR}/p${PREF}func_data.nii

                  python ${QCDIR}/save_image_diff.py -o ${OUTDIR} -i ${OUTDIR} -a ${PREF}func_data.nii -b p${PREF}func_data.nii -msg1 'Motion Corr' -msg2 'Motion & Physio Corr'

                  PREF=p${PREF}
                  mv ${OUTDIR}/*.png ${IMGDIR}
                  mv ${OUTDIR}/*.gif ${IMGDIR}
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

                  #if (( $(echo "$TR >= 2.5" | bc -l) )); then
                  #      # Used to preserve the full power spectrum
                  #      CF=1
                  #fi

                  # Correct slice acquisition to the middle of the volume acquisition
                  # To correct to the first slice, set REF_TIME to 0
                  # TODO [10.09.19] Pass this to the parameter section of the file
                  REF_TIME=`echo "${TR}/2" | bc -l`

                  #FIXME slice timing name is now not consistent
                  #print_debug 'Slice timing correction' "filtershift --in=${OUTDIR}/${PREF}func_data.nii \
                  #  --cf=${CF} --rt=${REF_TIME} --TR=${TR} --timing=${OUTDIR}/slice_timing_gs.txt --out=${OUTDIR}/a${PREF}func_data.nii"

                  #--rt=${REF_TIME}
                  # --rs=30 #--timing=${OUTDIR}/slice_timing_gs.txt
                  filtershift --in=${OUTDIR}/${PREF}func_data.nii --timing=${OUTDIR}/slice_acq.txt --cf=1 --timing=${OUTDIR}/slice_acq.txt --TR=${TR} --out=${OUTDIR}/a${PREF}func_data.nii
                  gunzip -f ${OUTDIR}/a${PREF}func_data.nii.gz
                  python ${QCDIR}/save_image_diff.py -o ${OUTDIR} -i ${OUTDIR} -a ${PREF}func_data.nii -b a${PREF}func_data.nii -msg1 'Before STC' -msg2 'After STC'
                  PREF=a${PREF}

                  mv ${OUTDIR}/*.png ${IMGDIR}
                  mv ${OUTDIR}/*.gif ${IMGDIR}
                  log_command_div 'Slice Timing correction'
            fi

            if [ "$DO_DPK" -eq "2" ]; then
                  print_debug 'Despiking [-localedit -NEW]'
                  3dDespike -nomask -NEW -localedit -prefix ${OUTDIR}/d${PREF}func_data.nii ${OUTDIR}/${PREF}func_data.nii
                  python ${QCDIR}/save_image_diff.py -o ${OUTDIR} -i ${OUTDIR} -a ${PREF}func_data.nii -b d${PREF}func_data.nii -msg1 'Before Despike' -msg2 'After Despike'
                  PREF=d${PREF}
                  mv ${OUTDIR}/*.png ${IMGDIR}
                  mv ${OUTDIR}/*.gif ${IMGDIR}
                  log_command_div 'Despike'
            fi



            # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
                        # X is odd
                        X=$((X+1))
                        HAS_ODD=1
                  fi

                  if (( $Y % 2 )); then
                        # Y is odd
                        Y=$((Y+1))
                        HAS_ODD=1
                  fi

                  #3dresample -prefix ${OUTDIR}/tmp.nii -orient ras -input ${OUTDIR}/rev_phase.nii
                  #3drefit -deoblique ${OUTDIR}/tmp.nii
                  #mv ${OUTDIR}/tmp.nii ${OUTDIR}/rev_phase.nii

                  if [ "${HAS_ODD}" -eq "1" ]; then
                        fslroi ${OUTDIR}/${PREF}func_data ${OUTDIR}/tmp 0 $X 0 $Y 0 $Z
                        gunzip -f ${OUTDIR}/tmp.nii.gz
                        mv ${OUTDIR}/tmp.nii  ${OUTDIR}/${PREF}func_data.nii

                        fslroi ${OUTDIR}/rev_phase ${OUTDIR}/tmp 0 $X 0 $Y 0 $Z
                        gunzip -f ${OUTDIR}/tmp.nii.gz
                        mv ${OUTDIR}/tmp.nii  ${OUTDIR}/rev_phase.nii

                  fi

                  if test -f "${OUTDIR}/voxel_displ_DZ.nii"; then
                        fslroi ${OUTDIR}/voxel_displ_DZ ${OUTDIR}/tmp 0 $X 0 $Y 0 $Z
                        gunzip -f ${OUTDIR}/tmp.nii.gz
                        mv ${OUTDIR}/tmp.nii  ${OUTDIR}/voxel_displ_DZ.nii

                        fslroi ${OUTDIR}/voxel_displ_DX ${OUTDIR}/tmp 0 $X 0 $Y 0 $Z
                        gunzip -f ${OUTDIR}/tmp.nii.gz
                        mv ${OUTDIR}/tmp.nii  ${OUTDIR}/voxel_displ_DX.nii

                        fslroi ${OUTDIR}/voxel_displ_DY ${OUTDIR}/tmp 0 $X 0 $Y 0 $Z
                        gunzip -f ${OUTDIR}/tmp.nii.gz
                        mv ${OUTDIR}/tmp.nii  ${OUTDIR}/voxel_displ_DY.nii
                  fi


                  # Use 5 volumes with normal and reverse phase encoding and create brain mask
                  HVOL=$((NVOLS/2))
                  VI=$((HVOL-2))
                  VF=$((HVOL+2))
                  3dTcat -tr ${TR} -prefix ${OUTDIR}/blipdown.nii ${OUTDIR}/${PREF}func_data.nii[${VI}..${VF}]
                  3dTcat -tr ${TR} -prefix ${OUTDIR}/blipup.nii ${OUTDIR}/rev_phase.nii

                  3dAutomask -apply_prefix ${OUTDIR}/blipup_m.nii ${OUTDIR}/blipup.nii
                  3dAutomask -apply_prefix ${OUTDIR}/blipdown_m.nii ${OUTDIR}/blipdown.nii

                  fslmaths ${OUTDIR}/blipup_m -nan -thrP 10 -ing 1000 ${OUTDIR}/blipup_mi
                  fslmaths ${OUTDIR}/blipdown_m -nan -thrP 10 -ing 1000 ${OUTDIR}/blipdown_mi

                  gunzip -f ${OUTDIR}/blipup_mi.nii.gz
                  gunzip -f ${OUTDIR}/blipdown_mi.nii.gz

                  copy_header ${OUTDIR}/blipup_mi.nii ${OUTDIR}/blipdown_mi.nii

                  #exit
                  # Concatenate both PE to create single 10-volume series from where the fieldmap will be esimated
                  3dTcat -tr ${TR} -prefix ${OUTDIR}/topup_data.nii ${OUTDIR}/blipdown_mi.nii[0..4] ${OUTDIR}/blipup_mi.nii[0..4]


                  rm -f ${OUTDIR}/blipdown.nii
                  rm -f ${OUTDIR}/blipup.nii
                  rm -f ${OUTDIR}/blipdown_m.nii
                  rm -f ${OUTDIR}/blipup_m.nii
                  rm -f ${OUTDIR}/blipdown_mi.nii
                  rm -f ${OUTDIR}/blipup_m.nii

                  python ./util/write_topup_file.py -out ${OUTDIR}/topupfield.txt

                  # Fieldmap estimation [ONLY estimation]
                  print_debug 'TOPUP command' "topup --imain=${OUTDIR}/topup_data.nii --datain=${OUTDIR}/topupfield.txt --out=${OUTDIR}/topup_results --fout=${OUTDIR}/fieldmap --iout=${OUTDIR}/w_topup_data --estmov=0,0,0  --regmod=membrane_energy  --minmet=1,1,1 --verbose --warpres=8,6,4 --miter=32,12,4 --subsamp=4,2,1 --fwhm=6,4,2"

                  #topup --imain=${OUTDIR}/topup_data.nii --datain=${OUTDIR}/topupfield.txt --out=${OUTDIR}/topup_results --fout=${OUTDIR}/fieldmap --iout=${OUTDIR}/w_topup_data \
                  #      --estmov=0,0,0  --regmod=membrane_energy  --minmet=1,1,1 --verbose --warpres=8,6,4 --miter=32,12,4 --subsamp=2,2,1 --fwhm=6,4,2

                  #applytopup --imain=${OUTDIR}/${PREF}func_data.nii --datain=${OUTDIR}/topupfield.txt --topup=${OUTDIR}/topup_results --inindex=1 --method=jac --interp=spline --out=${OUTDIR}/u${PREF}func_data
                  rm -f ${OUTDIR}/rtopup_data.nii
                  3dvolreg -prefix ${OUTDIR}/rtopup_data.nii ${OUTDIR}/topup_data.nii
                  copy_header ${OUTDIR}/topup_data.nii ${OUTDIR}/rtopup_data.nii
                  topup --imain=${OUTDIR}/rtopup_data.nii --datain=${OUTDIR}/topupfield.txt --out=${OUTDIR}/topup_results --fout=${OUTDIR}/fieldmap --iout=${OUTDIR}/w_topup_data \
                        --estmov=0,0,0  --regmod=bending_energy  --minmet=1,1,1 --verbose --warpres=16,12,8 --miter=16,10,5 --subsamp=2,2,1 --fwhm=16,12,8


                  applytopup --imain=${OUTDIR}/${PREF}func_data.nii --datain=${OUTDIR}/topupfield.txt --topup=${OUTDIR}/topup_results --inindex=1 --method=jac --interp=spline --out=${OUTDIR}/u${PREF}func_data
                  gunzip -f ${OUTDIR}/u${PREF}func_data.nii.gz

                  if ! test -f "${OUTDIR}/voxel_displ_DZ.nii"; then
                      applytopup --imain=${OUTDIR}/voxel_displ_DX.nii --datain=${OUTDIR}/topupfield.txt --topup=${OUTDIR}/topup_results --inindex=1 --method=jac --interp=spline --out=${OUTDIR}/wvoxel_displ_DX
                      mv ${OUTDIR}/wvoxel_displ_DX.nii ${OUTDIR}/voxel_displ_DX.nii

                      applytopup --imain=${OUTDIR}/voxel_displ_DY.nii --datain=${OUTDIR}/topupfield.txt --topup=${OUTDIR}/topup_results --inindex=1 --method=jac --interp=spline --out=${OUTDIR}/wvoxel_displ_DY
                      mv ${OUTDIR}/wvoxel_displ_DZ.nii ${OUTDIR}/voxel_displ_DZ.nii

                      applytopup --imain=${OUTDIR}/voxel_displ_DZ.nii --datain=${OUTDIR}/topupfield.txt --topup=${OUTDIR}/topup_results --inindex=1 --method=jac --interp=spline --out=${OUTDIR}/wvoxel_displ_DZ
                      mv ${OUTDIR}/wvoxel_displ_DZ.nii ${OUTDIR}/voxel_displ_DZ.nii
                  fi

                  python ${QCDIR}/save_image_diff.py -o ${OUTDIR} -i ${OUTDIR} -a ${PREF}func_data.nii -b u${PREF}func_data.nii -type mean -msg1 'Before FMAP' -msg2 'After FMAP'
                  PREF=u${PREF}
                  log_command_div 'FIELDMAP CORRECTION'
                  mv ${OUTDIR}/*.png ${IMGDIR}
                  mv ${OUTDIR}/*.gif ${IMGDIR}
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

                  if test -f "${OUTDIR}/voxel_displ_DX.nii"; then
                      fugue -i ${OUTDIR}/voxel_displ_DX.nii --dwell=${EES} --loadfmap=${OUTDIR}/fmap_phase_rads_func -u ${OUTDIR}/uvoxel_disp_DX.nii --unwarpdir=y-
                      fugue -i ${OUTDIR}/voxel_displ_DY.nii --dwell=${EES} --loadfmap=${OUTDIR}/fmap_phase_rads_func -u ${OUTDIR}/uvoxel_disp_DY.nii --unwarpdir=y-
                      fugue -i ${OUTDIR}/voxel_displ_DZ.nii --dwell=${EES} --loadfmap=${OUTDIR}/fmap_phase_rads_func -u ${OUTDIR}/uvoxel_disp_DZ.nii --unwarpdir=y-

                      mv ${OUTDIR}/uvoxel_disp_DX.nii ${OUTDIR}/voxel_displ_DX.nii
                      mv ${OUTDIR}/uvoxel_disp_DY.nii ${OUTDIR}/voxel_displ_DY.nii
                      mv ${OUTDIR}/uvoxel_disp_DZ.nii ${OUTDIR}/voxel_displ_DZ.nii
                  fi

                  gunzip -f ${OUTDIR}/u${PREF}func_data.nii.gz
                  python ${QCDIR}/save_image_diff.py -o ${OUTDIR} -i ${OUTDIR} -a ${PREF}func_data.nii -b u${PREF}func_data.nii -type mean -msg1 'Before FMAP' -msg2 'After FMAP'

                  PREF=u${PREF}
                  log_command_div 'FIELDMAP CORRECTION'
            fi

            # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


            if [ "$DO_DEOBL" -eq "1" ]; then
                  print_debug 'Deobliquing volumes'

                  #3dWarp -deoblique -newgrid ${DEOBL_VOX} -NN -prefix ${OUTDIR}/w${PREF}func_data.nii ${OUTDIR}/${PREF}func_data.nii
                  #cp ${OUTDIR}/${PREF}func_data.nii ${OUTDIR}/w${PREF}func_data.nii
                  3dresample -prefix ${OUTDIR}/w${PREF}func_data.nii -orient ras -input ${OUTDIR}/${PREF}func_data.nii
                  3drefit -deoblique ${OUTDIR}/w${PREF}func_data.nii


                  if ! test -f "${OUTDIR}/voxel_displ_DZ.nii"; then
                      3dresample -prefix ${OUTDIR}/wvoxel_displ_DX.nii -orient ras -input ${OUTDIR}/voxel_displ_DX.nii
                      3drefit -deoblique ${OUTDIR}/wvoxel_displ_DX.nii
                      mv ${OUTDIR}/wvoxel_displ_DX.nii ${OUTDIR}/voxel_displ_DX.nii

                      3dresample -prefix ${OUTDIR}/wvoxel_displ_DY.nii -orient ras -input ${OUTDIR}/voxel_displ_DY.nii
                      3drefit -deoblique ${OUTDIR}/wvoxel_displ_DY.nii
                      mv ${OUTDIR}/wvoxel_displ_DZ.nii ${OUTDIR}/voxel_displ_DZ.nii

                      3dresample -prefix ${OUTDIR}/wvoxel_displ_DZ.nii -orient ras -input ${OUTDIR}/voxel_displ_DZ.nii
                      3drefit -deoblique ${OUTDIR}/wvoxel_displ_DZ.nii
                      mv ${OUTDIR}/wvoxel_displ_DZ.nii ${OUTDIR}/voxel_displ_DZ.nii
                  fi

                  PREF=w${PREF}
                  log_command_div 'DEOBLIQUE'
            fi


            if [ "$DO_ING" -eq "1" ]; then
                  print_debug 'Global intensity normalization' 'fslmaths ${OUTDIR}/${PREF}func_data -ing ${GLOB_VAL} ${OUTDIR}/i${PREF}func_data'
                  fslmaths ${OUTDIR}/${PREF}func_data -nan -thrP 10 -ing ${GLOB_VAL} ${OUTDIR}/i${PREF}func_data
                  gunzip -f ${OUTDIR}/i${PREF}func_data.nii.gz
                  PREF=i${PREF}
            fi




            # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            # Update brain mask
            rm -f ${OUTDIR}/nat_mask.nii

            cp -f ${OUTDIR}/${PREF}func_data.nii ${OUTDIR}/proc_data_native.nii

            fslmaths ${OUTDIR}/proc_data_native.nii -Tmean ${OUTDIR}/proc_data_native_median.nii
            fast -n 4 --nopve -N -o ${OUTDIR}/func_seg.nii ${OUTDIR}/proc_data_native_median.nii
            fslmaths ${OUTDIR}/func_seg_seg.nii -thr 1.5 -bin -fillh26 -fillh26 -fillh -fillh ${OUTDIR}/nat_mask.nii
            gunzip -f ${OUTDIR}/nat_mask.nii

            rm -f ${OUTDIR}/proc_data_native_median.nii.gz


            #3dAutomask -prefix ${OUTDIR}/nat_mask.nii -apply_prefix ${OUTDIR}/proc_data_native.nii ${OUTDIR}/${PREF}func_data.nii
            #3dAutomask -prefix ${OUTDIR}/nat_mask.nii ${OUTDIR}/${PREF}func_data.nii


            3dcalc -a ${OUTDIR}/proc_data_native.nii -b ${OUTDIR}/nat_mask.nii -expr "a*b" -prefix ${OUTDIR}/proc_data_native_masked.nii

            mv ${OUTDIR}/proc_data_native_masked.nii ${OUTDIR}/proc_data_native.nii

            fslmaths ${OUTDIR}/proc_data_native.nii -Tmean ${OUTDIR}/mean_func_native.nii
            gunzip -f ${OUTDIR}/mean_func_native.nii.gz
            copy_header ${OUTDIR}/mean_func_native.nii ${OUTDIR}/nat_mask.nii


            3dcalc -a ${OUTDIR}/proc_data_native.nii -expr "a*1" -short -gscale -prefix ${OUTDIR}/proc_data_native_2.nii
            mv ${OUTDIR}/proc_data_native_2.nii  ${OUTDIR}/proc_data_native.nii

            if test -f "${OUTDIR}/voxel_displ_DZ.nii"; then
              3dcalc -a ${OUTDIR}/voxel_displ_DZ.nii -expr "a*1" -short -gscale -prefix ${OUTDIR}/voxel_displ_DZ_2.nii
              mv ${OUTDIR}/voxel_displ_DZ_2.nii  ${OUTDIR}/voxel_displ_DZ.nii

              3dcalc -a ${OUTDIR}/voxel_displ_DX.nii -expr "a*1" -short -gscale -prefix ${OUTDIR}/voxel_displ_DX_2.nii
              mv ${OUTDIR}/voxel_displ_DX_2.nii  ${OUTDIR}/voxel_displ_DX.nii

              3dcalc -a ${OUTDIR}/voxel_displ_DY.nii -expr "a*1" -short -gscale -prefix ${OUTDIR}/voxel_displ_DY_2.nii
              mv ${OUTDIR}/voxel_displ_DY_2.nii  ${OUTDIR}/voxel_displ_DY.nii

              copy_header ${OUTDIR}/proc_data_native.nii ${OUTDIR}/voxel_displ_DZ.nii
              copy_header ${OUTDIR}/proc_data_native.nii ${OUTDIR}/voxel_displ_DX.nii
              copy_header ${OUTDIR}/proc_data_native.nii ${OUTDIR}/voxel_displ_DY.nii
            fi


            mv ${OUTDIR}/*.png ${IMGDIR}
            mv ${OUTDIR}/*.gif ${IMGDIR}

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

                  #rm -f ${OUTDIR}/voxel_displ_DX.nii
                  #rm -f ${OUTDIR}/voxel_displ_DY.nii
                  #rm -f ${OUTDIR}/voxel_displ_DZ.nii
            fi


      } &> ${LOGDIR}/Native_PreProcessing.log
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

# Calculate WARP to group template
if [ "$DO_ANAT" -eq "1" ]; then
      printf "[$SUB] Functional <-> Anat <-> MNI Registration ... "
      START=$(date -u +%s.%N)


      {

          DOCOREG=1

          if [ "$DOCOREG" -eq "1" ]; then
                #ANATOMY processing
                ANATDIR=${OUTDIR}/anat
                mkdir ${ANATDIR}
                rm -f ${ANATDIR}/*

                cp ${OUTDIR}/t1.nii ${ANATDIR}/anat.nii

                3dresample -prefix ${ANATDIR}/tmp.nii -orient ras -input ${ANATDIR}/anat.nii
                mv ${ANATDIR}/tmp.nii ${ANATDIR}/anat.nii

                cp ${OUTDIR}/proc_data_native.nii ${ANATDIR}/func_data.nii

                NVOLS=`fslnvols ${ANATDIR}/func_data.nii`

                #FUNC
                fslmaths ${ANATDIR}/func_data.nii -Tmean -thrp 30 ${ANATDIR}/mean_func_data.nii
                gunzip -f ${ANATDIR}/mean_func_data.nii.gz
                ${ABIN}/DenoiseImage -i ${ANATDIR}/mean_func_data.nii -r 3x3x3 -n Gaussian -o [ ${ANATDIR}/mean_func_data_nds.nii ]


                #ANAT
                #3dresample -prefix ${ANATDIR}/mean_func_data_nds_.nii -orient ras -input ${ANATDIR}/mean_func_data_nds.nii
                #mv ${ANATDIR}/mean_func_data_nds_.nii ${ANATDIR}/mean_func_data_nds.nii
                #3drefit -deoblique ${ANATDIR}/mean_func_data_nds.nii

                #3dresample -prefix ${ANATDIR}/anat_orig.nii -orient ras -input ${ANATDIR}/anat.nii
                #3drefit -deoblique -duporigin ${ANATDIR}/mean_func_data_nds.nii ${ANATDIR}/anat_orig.nii

                ${ABIN}/N4BiasFieldCorrection -i ${ANATDIR}/anat.nii -s 3  -o [ ${ANATDIR}/anat_b.nii ]

                bet ${ANATDIR}/anat_b.nii ${ANATDIR}/anat_wb
                gunzip -f ${ANATDIR}/anat_wb.nii.gz

                #3dAutomask -apply_prefix ${ANATDIR}/anat_wb.nii ${ANATDIR}/anat_b.nii
                #3dUnifize -prefix ${ANATDIR}/anat_wb.nii  -input ${ANATDIR}/anat_b.nii
                #3dSkullStrip -input ${ANATDIR}/anat_wb.nii -prefix ${ANATDIR}/anat_wbb.nii

                #${ABIN}/ImageMath 3 ${ANATDIR}/anat_wb.nii RescaleImage ${ANATDIR}/anat_b.nii 0 1
                cp ${WDIR}/atlas/Prior/mni_icbm152_t1_brain_tal_nlin_sym_09c.nii ${ANATDIR}/template_t1.nii
                cp ${WDIR}/atlas/Prior/mni_icbm152_t1_tal_nlin_sym_09c_mask.nii ${ANATDIR}/brain_mask_mni.nii



                REF=${ANATDIR}/anat_wb.nii
                SRC=${ANATDIR}/template_t1.nii

                ${ABIN}/antsRegistration -d 3 -r [$REF,$SRC,1] -v 1 \
                            -m MI[$REF,$SRC,1,32] -t translation[0.1] -c [500,5.e-7,20] \
                            -s 3vox -f 3 -l 1 -n BSpline \
                            -m MI[$REF,$SRC,1,32,Regular,0.25] -t rigid[0.1] -c [500,5.e-7,20] \
                            -s 3vox -f 3 -l 1 -n BSpline \
                            -m MI[$REF,$SRC,1,32,Regular,0.25] -t affine[0.1] -c [500x100x10,5.e-7,10] \
                            -s 2x1x0vox -f 3x2x1 -l 1 -n BSpline \
                            -m CC[$REF,$SRC,1,3] -t SyN[0.2,3] -c [20,1.e-7,10] \
                            -s 4vox -f 3 -l 1 -n BSpline \
                            -o [${ANATDIR}/mni2anat,${ANATDIR}/mni2anat.nii]

                fslmaths ${ANATDIR}/mni2anat.nii -thr 10 -bin ${ANATDIR}/brain_mask.nii

                rm -f ${ANATDIR}/mni2anat1InverseWarp.nii.gz
                rm -f ${ANATDIR}/mni2anat1Warp.nii.gz
#                matlab "-nodesktop -nosplash " <<<"coreg_normalise('${ANATDIR}/template_t1.nii', {'${ANATDIR}/brain_mask_mni.nii'}, ${NVOLS}, '${ANATDIR}/anat_wb.nii', [2 0 0], {'${DARTEL_TEMPLATE_PREF}1.nii', '${DARTEL_TEMPLATE_PREF}2.nii', '${DARTEL_TEMPLATE_PREF}3.nii', '${DARTEL_TEMPLATE_PREF}4.nii', '${DARTEL_TEMPLATE_PREF}5.nii', '${DARTEL_TEMPLATE_PREF}6.nii'}, [5 5 5]); exit;"

                #python ${WDIR}/extract_brain.py -i "${ANATDIR}/anat_wb.nii" -o "${ANATDIR}/brain_mask.nii"
                #exit
                ${ABIN}/DenoiseImage -i ${ANATDIR}/anat_b.nii -p 1x1x1 -r 4x4x4 -n Rician -x ${ANATDIR}/brain_mask.nii -o [ ${ANATDIR}/anat_wbn.nii ]
                3dcalc -a ${ANATDIR}/anat_wbn.nii -b ${ANATDIR}/brain_mask.nii -expr 'a*b' -prefix ${ANATDIR}/anat_proc_brain.nii


                mv ${ANATDIR}/anat_wbn.nii ${ANATDIR}/anat_proc.nii
                copy_header ${ANATDIR}/anat.nii ${ANATDIR}/anat_proc.nii

                rm -f ${ANATDIR}/anat_w*.nii
                rm -f ${ANATDIR}/anat.nii
                rm -f ${ANATDIR}/template_t1.nii

                REF=${ANATDIR}/anat_proc.nii
                SRC=${ANATDIR}/mean_func_data_nds.nii

                ${ABIN}/antsRegistration -d 3 -r [$REF,$SRC,1] -v 1 \
                            -m MI[$REF,$SRC,1,32] -t translation[0.1] -c [500,5.e-7,20] \
                            -s 3vox -f 3 -l 1 -n BSpline \
                            -o [${ANATDIR}/anat_proc_a,${ANATDIR}/anat_proc_a.nii]

                ${ABIN}/antsApplyTransforms \
                      -i ${ANATDIR}/anat_proc.nii \
                      --float \
                      -r ${ANATDIR}/anat_proc.nii \
                      -o ${ANATDIR}/anat_proc_a2.nii \
                      -t [${ANATDIR}/anat_proc_a0GenericAffine.mat,1] \
                      -n NearestNeighbor

                 ${ABIN}/antsApplyTransforms \
                       -i ${ANATDIR}/anat_proc_brain.nii \
                       --float \
                       -r ${ANATDIR}/anat_proc.nii \
                       -o ${ANATDIR}/anat_proc_a3.nii \
                       -t [${ANATDIR}/anat_proc_a0GenericAffine.mat,1] \
                       -n NearestNeighbor

                mv ${ANATDIR}/anat_proc_a2.nii ${ANATDIR}/anat_proc.nii
                mv ${ANATDIR}/anat_proc_a3.nii ${ANATDIR}/anat_proc_brain.nii
                rm -f ${ANATDIR}/anat_proc_a.nii



                #copy_header ${ANATDIR}/mean_func_data.nii ${ANATDIR}/mean_func_data_nds.nii

                matlab "-nodesktop -nosplash " <<<"coreg_normalise('${ANATDIR}/mean_func_data_nds.nii', '${OUTDIR}/proc_data_native.nii', ${NVOLS}, '${ANATDIR}/anat_proc_brain.nii', [1 0 0], {'${DARTEL_TEMPLATE_PREF}1.nii', '${DARTEL_TEMPLATE_PREF}2.nii', '${DARTEL_TEMPLATE_PREF}3.nii', '${DARTEL_TEMPLATE_PREF}4.nii', '${DARTEL_TEMPLATE_PREF}5.nii', '${DARTEL_TEMPLATE_PREF}6.nii'}, [5 5 5]); exit;"

                copy_header ${ANATDIR}/anat_proc_brain.nii ${ANATDIR}/anat_proc.nii

                matlab "-nodesktop -nosplash " <<<"coreg_normalise('${ANATDIR}/mean_func_data_nds.nii', '${OUTDIR}/proc_data_native.nii', ${NVOLS}, '${ANATDIR}/anat_proc.nii', [0 1 0], {'${DARTEL_TEMPLATE_PREF}1.nii', '${DARTEL_TEMPLATE_PREF}2.nii', '${DARTEL_TEMPLATE_PREF}3.nii', '${DARTEL_TEMPLATE_PREF}4.nii', '${DARTEL_TEMPLATE_PREF}5.nii', '${DARTEL_TEMPLATE_PREF}6.nii'}, [5 5 5]); exit;"

                #matlab "-nodesktop -nosplash " <<<"affine_spm2fsl('${ANATDIR}', '${ANATDIR}/anat_proc.nii','${ANATDIR}/mean_func_data_nds.nii'); exit;"

          fi
          rm -f ${ANATDIR}/anat_native.nii

          mv ${ANATDIR}/rc1anat_proc.nii ${ANATDIR}/rc1anat_proc_dartel.nii
          mv ${ANATDIR}/rc2anat_proc.nii ${ANATDIR}/rc2anat_proc_dartel.nii
          mv ${ANATDIR}/rc3anat_proc.nii ${ANATDIR}/rc3anat_proc_dartel.nii

          matlab "-nodesktop -nosplash " <<<"reslice_images('${ANATDIR}/mean_func_data_nds.nii', {'${ANATDIR}/anat_proc.nii', '${ANATDIR}/anat_proc_brain.nii', '${ANATDIR}/c1anat_proc.nii', '${ANATDIR}/c2anat_proc.nii', '${ANATDIR}/c3anat_proc.nii'}); exit;"
          #exit

          mv ${ANATDIR}/ranat_proc.nii ${ANATDIR}/anat_native.nii
          mv ${ANATDIR}/ranat_proc_brain.nii ${ANATDIR}/anat_native_brain.nii

          mv ${ANATDIR}/rc1anat_proc.nii ${ANATDIR}/gm_tpm_native.nii
          mv ${ANATDIR}/rc2anat_proc.nii ${ANATDIR}/wm_tpm_native.nii
          mv ${ANATDIR}/rc3anat_proc.nii ${ANATDIR}/csf_tpm_native.nii


          3dTcat -prefix ${OUTDIR}/example_func.nii ${OUTDIR}/proc_data_native.nii[0]

          copy_header ${ANATDIR}/mean_func_data_nds.nii ${ANATDIR}/anat_native.nii
          copy_header ${ANATDIR}/mean_func_data_nds.nii ${ANATDIR}/anat_native_brain.nii
          copy_header ${ANATDIR}/mean_func_data_nds.nii ${ANATDIR}/gm_tpm_native.nii
          copy_header ${ANATDIR}/mean_func_data_nds.nii ${ANATDIR}/wm_tpm_native.nii
          copy_header ${ANATDIR}/mean_func_data_nds.nii ${ANATDIR}/csf_tpm_native.nii



          flirt -in ${ANATDIR}/anat_native.nii -ref ${ANATDIR}/anat_proc.nii -omat ${ANATDIR}/func2anat_fsl.mat -out ${ANATDIR}/anat_native_fsl \
               -usesqform -dof 6 -coarsesearch 18 -finesearch 9 -searchrx -20 20 -searchry -20 20 -searchrz -20 20


          convert_xfm -omat ${ANATDIR}/anat2func_fsl.mat -inverse ${ANATDIR}/func2anat_fsl.mat
          flirt -applyxfm -init ${ANATDIR}/anat2func_fsl.mat -in ${ANATDIR}/anat_proc.nii -out ${ANATDIR}/anat_native_fsl_func.nii.gz -ref ${ANATDIR}/mean_func_data_nds.nii

          gunzip -f ${ANATDIR}/anat_native_fsl_func.nii.gz
          copy_header ${ANATDIR}/mean_func_data_nds.nii ${ANATDIR}/anat_native_fsl_func.nii
          gzip ${ANATDIR}/anat_native_fsl_func.nii
          rm -f ${ANATDIR}/anat_native_fsl_func.nii

          fslmaths ${ANATDIR}/wm_tpm_native.nii -add ${ANATDIR}/gm_tpm_native.nii -thr 0.1 -bin -fillh -nan ${OUTDIR}/nat_mask_tpm.nii
          gunzip -f ${OUTDIR}/nat_mask_tpm.nii.gz
      if [ "$DO_ANAT_FMAP" -eq "1"]; then

            REF=${ANATDIR}/anat_native_brain.nii
            SRC=${ANATDIR}/mean_func_data_nds.nii
            ${ABIN}/antsRegistration -d 3 -r [$REF,$SRC,0] -v 1 \
                        -m CC[$REF,$SRC,1,3] -t SyN[0.05,1] -c [20x10,1.e-7,10] \
                        -s 6x4vox -f 4x2 -l 1 -n BSpline  -g 0x1x0 \
                        -o [${ANATDIR}/fmap,${ANATDIR}/u_mean_func_data_nds.nii]

            copy_header ${ANATDIR}/mean_func_data_nds.nii  ${ANATDIR}/u_mean_func_data_nds.nii
            NVOLS=`fslnvols ${OUTDIR}/proc_data_native.nii`
            mkdir ${ANATDIR}/tmp
            3dTsplit4D -prefix ${ANATDIR}/tmp/func_data.nii -keep_datum ${OUTDIR}/proc_data_native.nii

            parallel -j4 --line-buffer correct_anat_fmap ::: $(seq 0 ${NVOLS}) ::: ${ANATDIR} ::: ${ABIN} ::: fmap

            rm -f ${ANATDIR}/proc_data_native_u.nii
            3dTcat -prefix ${OUTDIR}/proc_data_native_u.nii -tr ${TR} ${ANATDIR}/tmp/func_data_n*.nii

            copy_header ${OUTDIR}/proc_data_native.nii  ${OUTDIR}/proc_data_native_u.nii

            rm -f -r ${ANATDIR}/tmp

      fi
      if [ "$DO_CLEAN" -eq "1" ]; then

            rm -f ${ANATDIR}/c1anat_proc.nii
            rm -f ${ANATDIR}/c2anat_proc.nii
            rm -f ${ANATDIR}/c3anat_proc.nii
            rm -f ${ANATDIR}/c4anat_proc.nii
            rm -f ${ANATDIR}/c5anat_proc.nii
            rm -f ${ANATDIR}/c6anat_proc.nii

            rm -f ${ANATDIR}/anat.nii
            rm -f ${ANATDIR}/anat_b.nii
            rm -f ${ANATDIR}/func_data.nii
            rm -f ${ANATDIR}/func_data.mat
            rm -f ${ANATDIR}/anat_wb.nii
            rm -f ${ANATDIR}/brain_mask.nii
            gzip -f ${ANATDIR}/iy_anat_proc.nii
            gzip -f ${ANATDIR}/y_anat_proc.nii

      fi
    } &> ${LOGDIR}/T1_Processing.log

      END=$(date -u +%s.%N)
      DIFF=`echo "( $END - $START )" | bc`
      printf "DONE [%.1f s]\n" $DIFF
fi






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

      rm -f -r ${OUTDIR}/tmp/

      DOIC=${RUN_ICA}
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

            rm -f ${OUTDIR}/tmp_melodic.nii
            # Create a temporary detrended 3d+time series to run melodic
            # For ICA-FIX, smoothing is not recommended [only during visualization]
            rm -f ${OUTDIR}/tmp_melodic.nii
            rm -f ${OUTDIR}/nat_mask.nii
            #if test -f ${OUTDIR}/phyca_data__NN_map.nii; then
            #  #cp ${OUTDIR}/phyca_data__NN_map.nii ${OUTDIR}/nat_mask_dil.nii
            #  #3dcalc -a ${OUTDIR}/proc_data_native.nii -expr "astep(a,0)" -prefix ${OUTDIR}/nat_mask_dil.nii
            #  fslmaths ${OUTDIR}/proc_data_native.nii -Tmean -bin ${OUTDIR}/nat_mask_dil.nii
            #  gunzip -f ${OUTDIR}/nat_mask_dil.nii.gz
            #else
            #3dAutomask -prefix ${OUTDIR}/nat_mask.nii -peel 2  ${OUTDIR}/proc_data_native.nii
            #fi


            copy_header ${ANATDIR}/mean_func_data_nds.nii ${OUTDIR}/nat_mask_tpm.nii
            # -stopband 0 0.009
            3dTproject -prefix ${OUTDIR}/tmp_melodic.nii -polort 2  -mask ${OUTDIR}/nat_mask_tpm.nii -TR ${TR} -input ${OUTDIR}/proc_data_native.nii

            rm -f -R ${OUTDIR}/melodic.ic
            rm -f -R ${OUTDIR}/FIX

            copy_header ${OUTDIR}/proc_data_native.nii ${OUTDIR}/tmp_melodic.nii
            #-vnorm
            3dpc -eigonly -nscale -vmean -vnorm -mask ${OUTDIR}/nat_mask_tpm.nii -prefix ${OUTDIR}/pc_var ${OUTDIR}/tmp_melodic.nii

            NIC=`python ${WDIR}/util/ica_model_order.py -outfile ${OUTDIR}/num_ics.txt -eig_file ${OUTDIR}/pc_var_eig.1D`

            #NIC=`tail -n 1 ${OUTDIR}/num_ics.txt`

            # --dimest=lap is the default, but it seems empirically to overestimate components
            # see also Varoquaux et al. 2010 NeuroImage
            if [ "$ICA_TYPE" == "melodic" ]; then
                  melodic -i ${OUTDIR}/tmp_melodic.nii -o ${OUTDIR}/melodic.ic --tr=${TR}  --mmthresh=0.5 --nobet --dim=${NIC}  \
                         --mask=${OUTDIR}/nat_mask_tpm.nii --nl=pow4 --Ostats --report --eps=0.0001 --bgimage=${ANATDIR}/anat_native_fsl_func
            fi

            if [ "$ICA_TYPE" == "canica" ]; then
                  mkdir ${OUTDIR}/melodic.ic
                  #mkdir ${OUTDIR}/melodic.ic/report

                  #if test -f "${OUTDIR}/Group_Comps_Native.nii"; then
                  #      python ${WDIR}/canica_estimation.py -o ${OUTDIR}/melodic.ic -in ${OUTDIR}/tmp_melodic.nii \
                  #              -mask ${OUTDIR}/nat_mask.nii -nIC ${NIC} -decomp dictlearning -dict_init ${OUTDIR}/Group_Comps_Native.nii
                  #else
                        python ${WDIR}/canica_estimation.py -o ${OUTDIR}/melodic.ic -in ${OUTDIR}/tmp_melodic.nii \
                                -mask ${OUTDIR}/nat_mask_tpm.nii -nIC ${NIC} #-decomp dictlearning
                  #fi

                  copy_header ${OUTDIR}/nat_mask.nii ${OUTDIR}/melodic.ic/CanICA_Comps_Z.nii
                  copy_header ${OUTDIR}/nat_mask.nii ${OUTDIR}/melodic.ic/CanICA_Comps.nii

                  mv ${OUTDIR}/melodic.ic/CanICA_Comps_Z.nii ${OUTDIR}/melodic.ic/melodic_IC.nii
                  gzip ${OUTDIR}/melodic.ic/melodic_IC.nii
                  rm ${OUTDIR}/melodic.ic/melodic_IC.nii

                  # Remove NaNs from the IC maps, in case they are there
                  fslmaths ${OUTDIR}/melodic.ic/melodic_IC.nii.gz -nan -mas ${OUTDIR}/nat_mask.nii ${OUTDIR}/melodic.ic/melodic_IC.nii.gz

                  gunzip -f ${OUTDIR}/melodic.ic/melodic_IC.nii.gz
                  copy_header ${OUTDIR}/melodic.ic/CanICA_Comps.nii ${OUTDIR}/melodic.ic/melodic_IC.nii
                  gzip -f ${OUTDIR}/melodic.ic/melodic_IC.nii

                  echo "1" > ${OUTDIR}/melodic.ic/grot.txt
                  melodic -i ${OUTDIR}/melodic.ic/melodic_IC.nii.gz --ICs=${OUTDIR}/melodic.ic/melodic_IC.nii --nobet \
                      --nomask --bgthreshold=0.001 --tr=${TR}   \
                      --mix=${OUTDIR}/melodic.ic/grot.txt -o ${OUTDIR}/melodic.ic --Oall --report -v --mmthresh=0.5

                  # usage: extract_ics.py [-h] -o OUTDIR -i INFILE -icmap ICMAP
                  python ${WDIR}/util/extract_ics.py -i ${OUTDIR}/proc_data_native.nii  \
                          -icmap ${OUTDIR}/melodic.ic/CanICA_Comps.nii -nic ${NIC} -out ${OUTDIR}/melodic.ic/  #-decomp dictlearning

            fi


            rm ${OUTDIR}/tmp_melodic.nii

            # Creates a smoothed copy of the melodic components. For visualization purposes only
            3dmerge -doall -1blur_fwhm 5 -prefix ${OUTDIR}/melodic.ic/melodic_IC_smooth.nii ${OUTDIR}/melodic.ic/melodic_IC.nii.gz
            gunzip -f ${OUTDIR}/melodic.ic/melodic_IC.nii.gz
            copy_header ${OUTDIR}/melodic.ic/melodic_IC.nii ${OUTDIR}/melodic.ic/melodic_IC_smooth.nii
            gzip ${OUTDIR}/melodic.ic/melodic_IC.nii


      fi


            mkdir ${OUTDIR}/tmp
            #cp -f ${FIXDIR}/hand_labels_noise.txt ${OUTDIR}/hand_labels_noise.txt
            cp -f -r ${FIXDIR}/filtered_func_data.ica/* ${OUTDIR}/tmp
            rm -f -r ${FIXDIR}
            mkdir ${FIXDIR}
            mkdir ${FIXDIR}/mc
            mkdir ${FIXDIR}/reg
            mkdir ${FIXDIR}/filtered_func_data.ica

            cp -f -r ${OUTDIR}/tmp/* ${FIXDIR}/filtered_func_data.ica
            rm -f -r ${OUTDIR}/tmp

            cp -f ${ANATDIR}/anat_proc_brain.nii ${FIXDIR}/reg/highres.nii
            gzip -f ${FIXDIR}/reg/highres.nii

            cp -f ${OUTDIR}/proc_data_native.nii ${FIXDIR}/filtered_func_data.nii
            gzip -f ${FIXDIR}/filtered_func_data.nii
            cp -f ${OUTDIR}/motion_estimate_fsl.par ${FIXDIR}/mc/prefiltered_func_data_mcf.par
            cp -f -R ${OUTDIR}/melodic.ic/* ${FIXDIR}/filtered_func_data.ica
            cp ${OUTDIR}/nat_mask_tpm.nii ${FIXDIR}/mask.nii
            gzip -f ${FIXDIR}/mask.nii


            fslmaths ${OUTDIR}/proc_data_native.nii -Tmean  ${FIXDIR}/reg/example_func.nii
            cp -f ${FIXDIR}/reg/example_func.nii.gz ${FIXDIR}/mean_func.nii.gz

            #cp -f ${ANATDIR}/mean_func_data_nds.nii ${FIXDIR}/reg/example_func.nii
            #cp -f ${ANATDIR}/mean_func_data_nds.nii ${FIXDIR}/mean_func.nii
            #gzip -f ${FIXDIR}/reg/example_func.nii
            #gzip -f ${FIXDIR}/mean_func.nii

            rm -f  ${FIXDIR}/reg/highres2example_func.mat
            cp -f ${ANATDIR}/anat2func_fsl.mat ${FIXDIR}/reg/highres2example_func.mat

            gunzip -f ${FIXDIR}/mean_func.nii.gz
            gunzip -f ${FIXDIR}/filtered_func_data.ica/melodic_IC.nii.gz

            copy_header ${FIXDIR}/mean_func.nii ${FIXDIR}/filtered_func_data.ica/melodic_IC_smooth.nii
            copy_header ${FIXDIR}/mean_func.nii ${FIXDIR}/filtered_func_data.ica/melodic_IC.nii

            gzip -f ${FIXDIR}/mean_func.nii
            gzip -f ${FIXDIR}/filtered_func_data.ica/melodic_IC.nii

            rm -f ${FIXDIR}/reg/highres2exfunc.nii.gz
            mv ${FIXDIR}/hand_labels_noise.txt ${OUTDIR}/hand_labels_noise.txt


            # HERE RUN THIS AFTERWARDS
            CURDIR=`pwd`
            cd ${FIXBIN}
            source ${FIXBIN}/fix -f ${FIXDIR}
            cd ${CURDIR}


            gzip -f ${FIXDIR}/filtered_func_data.nii
            gzip -f ${FIXDIR}/mask.nii
            rm -r -f  ${OUTDIR}/melodic.ic
      } &> ${LOGDIR}/ICA.log

      END=$(date -u +%s.%N)
      DIFF=`echo "( $END - $START )" | bc`
      printf "DONE [%.1f s]\n" $DIFF

fi

# Apply normalisation WARP
# TODO --> Change Name to something more meaningful
if [ "$DO_MNI" -eq "1" ]; then

      printf "[$SUB] Converting Atlas to Native Space ... "
      START=$(date -u +%s.%N)

      {
            WDIR=$(awk -F\=  '/^WDIR/{print $2}' "${CONFIG}" | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
            cd "${WDIR}"

            rm -f ${ANATDIR}/c6anat_proc.nii


            rm -f ${ANATDIR}/lg400_cobra_group.nii
            rm -f ${ANATDIR}/yeo17cobra_group.nii
            rm -f ${ANATDIR}/template_group.nii

#NII=(cingulo_opercular default_mode dorsal_attention dorsal_somatomotor \
#      ventral_attention early_auditory language lateral_prefrontal \
#      left_frontoparietal medial_prefrontal right_frontoparietal \
#      ventral_somatomotor visual_parafoveal visual_peripheral)
#
#BASE=/home/luna.kuleuven.be/u0101486/workspace/data/ConnectEx/Template/
#3dTcat -prefix ${BASE}/mantini2013.nii ${BASE}/${NII[0]}_group.nii ${BASE}/${NII[1]}_group.nii \
#      ${BASE}/${NII[2]}_group.nii ${BASE}/${NII[3]}_group.nii ${BASE}/${NII[4]}_group.nii \
#      ${BASE}/${NII[5]}_group.nii ${BASE}/${NII[6]}_group.nii ${BASE}/${NII[7]}_group.nii \
#      ${BASE}/${NII[8]}_group.nii ${BASE}/${NII[9]}_group.nii ${BASE}/${NII[10]}_group.nii \
#      ${BASE}/${NII[11]}_group.nii ${BASE}/${NII[12]}_group.nii ${BASE}/${NII[13]}_group.nii
#
#exit
#NII_IN=(    \
#         \
#        )


            # TODO Create a 4D map for the Mantini Maps
            cp ${GROUPDIR}/lg400_cobra_group.nii ${ANATDIR}/lg400_cobra_group.nii
            cp ${GROUPDIR}/yeo17cobra_group.nii ${ANATDIR}/yeo17cobra_group.nii

            cp ${GROUPDIR}/Group_T1_Avg_Brain.nii ${ANATDIR}/template_group.nii

            if test -f ${GROUPDIR}/cingulo_opercular_group.nii; then
              cp ${GROUPDIR}/cingulo_opercular_group.nii ${ANATDIR}/cingulo_opercular_group.nii
              cp ${GROUPDIR}/default_mode_group.nii ${ANATDIR}/default_mode_group.nii
              cp ${GROUPDIR}/dorsal_attention_group.nii ${ANATDIR}/dorsal_attention_group.nii
              cp ${GROUPDIR}/dorsal_somatomotor_group.nii ${ANATDIR}/dorsal_somatomotor_group.nii
              cp ${GROUPDIR}/ventral_attention_group.nii ${ANATDIR}/ventral_attention_group.nii
              cp ${GROUPDIR}/early_auditory_group.nii ${ANATDIR}/early_auditory_group.nii
              cp ${GROUPDIR}/language_group.nii ${ANATDIR}/language_group.nii
              cp ${GROUPDIR}/lateral_prefrontal_group.nii ${ANATDIR}/lateral_prefrontal_group.nii
              cp ${GROUPDIR}/left_frontoparietal_group.nii ${ANATDIR}/left_frontoparietal_group.nii
              cp ${GROUPDIR}/medial_prefrontal_group.nii ${ANATDIR}/medial_prefrontal_group.nii
              cp ${GROUPDIR}/right_frontoparietal_group.nii ${ANATDIR}/right_frontoparietal_group.nii
              cp ${GROUPDIR}/ventral_somatomotor_group.nii ${ANATDIR}/ventral_somatomotor_group.nii
              cp ${GROUPDIR}/visual_parafoveal_group.nii ${ANATDIR}/visual_parafoveal_group.nii
              cp ${GROUPDIR}/visual_peripheral_group.nii ${ANATDIR}/visual_peripheral_group.nii
            fi

            if ! test -f ${ANATDIR}/u_mean_func_data_nds.nii; then
              cp ${ANATDIR}/mean_func_data_nds.nii ${ANATDIR}/u_mean_func_data_nds.nii
            fi

            if test -f ${GROUPDIR}/cingulo_opercular_group.nii; then
              MANTININETS="'${ANATDIR}/cingulo_opercular_group.nii',"
              MANTININETS=${MANTININETS}"'${ANATDIR}/default_mode_group.nii',"
              MANTININETS=${MANTININETS}"'${ANATDIR}/dorsal_attention_group.nii',"
              MANTININETS=${MANTININETS}"'${ANATDIR}/dorsal_somatomotor_group.nii',"
              MANTININETS=${MANTININETS}"'${ANATDIR}/ventral_attention_group.nii',"
              MANTININETS=${MANTININETS}"'${ANATDIR}/early_auditory_group.nii',"
              MANTININETS=${MANTININETS}"'${ANATDIR}/language_group.nii',"
              MANTININETS=${MANTININETS}"'${ANATDIR}/lateral_prefrontal_group.nii',"
              MANTININETS=${MANTININETS}"'${ANATDIR}/left_frontoparietal_group.nii',"
              MANTININETS=${MANTININETS}"'${ANATDIR}/medial_prefrontal_group.nii',"
              MANTININETS=${MANTININETS}"'${ANATDIR}/right_frontoparietal_group.nii',"
              MANTININETS=${MANTININETS}"'${ANATDIR}/ventral_somatomotor_group.nii',"
              MANTININETS=${MANTININETS}"'${ANATDIR}/visual_parafoveal_group.nii',"
              MANTININETS=${MANTININETS}"'${ANATDIR}/visual_peripheral_group.nii'"
              matlab "-nodesktop -nosplash " <<<"coreg_normalise('${ANATDIR}/u_mean_func_data_nds.nii', {'${ANATDIR}/lg400_cobra_group.nii', '${ANATDIR}/yeo17cobra_group.nii', ${MANTININETS}}, 1, '${ANATDIR}/anat_proc.nii', [0 0 2], {'${DARTEL_TEMPLATE_PREF}1.nii', '${DARTEL_TEMPLATE_PREF}2.nii', '${DARTEL_TEMPLATE_PREF}3.nii', '${DARTEL_TEMPLATE_PREF}4.nii', '${DARTEL_TEMPLATE_PREF}5.nii', '${DARTEL_TEMPLATE_PREF}6.nii'}, [0 0 0], 0); exit;"
            else
              matlab "-nodesktop -nosplash " <<<"coreg_normalise('${ANATDIR}/u_mean_func_data_nds.nii', {'${ANATDIR}/lg400_cobra_group.nii', '${ANATDIR}/yeo17cobra_group.nii'}, 1, '${ANATDIR}/anat_proc.nii', [0 0 2], {'${DARTEL_TEMPLATE_PREF}1.nii', '${DARTEL_TEMPLATE_PREF}2.nii', '${DARTEL_TEMPLATE_PREF}3.nii', '${DARTEL_TEMPLATE_PREF}4.nii', '${DARTEL_TEMPLATE_PREF}5.nii', '${DARTEL_TEMPLATE_PREF}6.nii'}, [0 0 0], 0); exit;"
            fi
            #matlab "-nodesktop -nosplash " <<<"coreg_normalise('${ANATDIR}/u_mean_func_data_nds.nii', {'${ANATDIR}/mantini2013_group.nii'}, 14, '${ANATDIR}/anat_proc.nii', [0 0 2], {'${DARTEL_TEMPLATE_PREF}1.nii', '${DARTEL_TEMPLATE_PREF}2.nii', '${DARTEL_TEMPLATE_PREF}3.nii', '${DARTEL_TEMPLATE_PREF}4.nii', '${DARTEL_TEMPLATE_PREF}5.nii', '${DARTEL_TEMPLATE_PREF}6.nii'}, [0 0 0], 0); exit;"


            matlab "-nodesktop -nosplash " <<<"coreg_same_image('${ANATDIR}/u_mean_func_data_nds.nii', '${ANATDIR}/wlg400_cobra_group_u_rc1anat_proc_dartel.nii', '', 0); exit;"
            matlab "-nodesktop -nosplash " <<<"coreg_same_image('${ANATDIR}/u_mean_func_data_nds.nii', '${ANATDIR}/wyeo17cobra_group_u_rc1anat_proc_dartel.nii', '', 0); exit;"

            if test -f ${GROUPDIR}/cingulo_opercular_group.nii; then

              NII_IN=(cingulo_opercular default_mode dorsal_attention dorsal_somatomotor \
                    ventral_attention early_auditory language lateral_prefrontal \
                    left_frontoparietal medial_prefrontal right_frontoparietal \
                    ventral_somatomotor visual_parafoveal visual_peripheral)

              for NII in "${NII_IN[@]}";
              do
              #  cp ${WDIR}/atlas/Mantini2013/RSNs_fMRI/${NII}.nii ${GROUPDIR}/${NII}_mni.nii
                 matlab "-nodesktop -nosplash " <<<"coreg_same_image('${ANATDIR}/u_mean_func_data_nds.nii', '${ANATDIR}/w${NII}_group_u_rc1anat_proc_dartel.nii', '', 0); exit;"
                 copy_header ${ANATDIR}/mean_func_data_nds.nii  ${ANATDIR}/rw${NII}_group_u_rc1anat_proc_dartel.nii
                 mv ${ANATDIR}/rw${NII}_group_u_rc1anat_proc_dartel.nii ${ANATDIR}/${NII}_native.nii
              done

              NII_IN=(cingulo_opercular default_mode dorsal_attention dorsal_somatomotor \
                    ventral_attention early_auditory language lateral_prefrontal \
                    left_frontoparietal medial_prefrontal right_frontoparietal \
                    ventral_somatomotor visual_parafoveal visual_peripheral)

              3dTcat -prefix ${ANATDIR}/mantini2013_native.nii ${ANATDIR}/${NII_IN[0]}_native.nii \
                  ${ANATDIR}/${NII_IN[1]}_native.nii ${ANATDIR}/${NII_IN[2]}_native.nii ${ANATDIR}/${NII_IN[3]}_native.nii \
                  ${ANATDIR}/${NII_IN[4]}_native.nii ${ANATDIR}/${NII_IN[5]}_native.nii ${ANATDIR}/${NII_IN[6]}_native.nii \
                  ${ANATDIR}/${NII_IN[7]}_native.nii ${ANATDIR}/${NII_IN[8]}_native.nii ${ANATDIR}/${NII_IN[9]}_native.nii \
                  ${ANATDIR}/${NII_IN[10]}_native.nii ${ANATDIR}/${NII_IN[11]}_native.nii ${ANATDIR}/${NII_IN[12]}_native.nii \
                  ${ANATDIR}/${NII_IN[13]}_native.nii

              for NII in "${NII_IN[@]}";
              do
                rm -f ${ANATDIR}/${NII}_native.nii
                rm -f ${ANATDIR}/${NII}_group.nii
                rm -f ${ANATDIR}/w${NII}_group_u_rc1anat_proc_dartel.nii
              done

            fi


            copy_header ${ANATDIR}/mean_func_data_nds.nii  ${ANATDIR}/rwlg400_cobra_group_u_rc1anat_proc_dartel.nii
            copy_header ${ANATDIR}/mean_func_data_nds.nii  ${ANATDIR}/rwyeo17cobra_group_u_rc1anat_proc_dartel.nii

            mv ${ANATDIR}/rwlg400_cobra_group_u_rc1anat_proc_dartel.nii ${ANATDIR}/lg400_cobra_native.nii
            mv ${ANATDIR}/rwyeo17cobra_group_u_rc1anat_proc_dartel.nii ${ANATDIR}/yeo17cobra_native.nii

            rm -f ${ANATDIR}/wlg400_cobra_group_u_rc1anat_proc_dartel.nii
            rm -f ${ANATDIR}/wyeo17cobra_group_u_rc1anat_proc_dartel.nii

            rm -f ${ANATDIR}/lg400_cobra_group.nii
            rm -f ${ANATDIR}/yeo17cobra_group.nii
            rm -f ${ANATDIR}/template_group.nii



} &> ${LOGDIR}/MNI_to_Native.log

      END=$(date -u +%s.%N)
      DIFF=`echo "( $END - $START )" | bc`
      printf "DONE [%.1f s]\n" $DIFF
fi


# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------


if  [ "${EXTRACT_NUIS}" -eq "1" ]; then
      printf "[$SUB] Extracting Nuisance Regressors ... "
      START=$(date -u +%s.%N)
      {

      SUFFIX=''
      if test -f "${OUTDIR}/proc_data_native_nlm.nii"; then
            SUFFIX='_nlm'
      fi
      #copy_header ${ANATDIR}/mean_func_nds.nii  ${ANATDIR}/csf_tpm_native.nii
      #copy_header ${ANATDIR}/mean_func_nds.nii  ${ANATDIR}/gm_tpm_native.nii
      #copy_header ${ANATDIR}/mean_func_nds.nii  ${ANATDIR}/wm_tpm_native.nii

      # ERODE tissue masks
      fslmaths ${ANATDIR}/csf_tpm_native.nii -thr 0.75 -bin ${ANATDIR}/csf_mask_native.nii
      gunzip -f ${ANATDIR}/csf_mask_native.nii.gz

      fslmaths ${ANATDIR}/wm_tpm_native.nii -thr 0.75 -bin ${ANATDIR}/wm_mask_native.nii
      gunzip -f ${ANATDIR}/wm_mask_native.nii.gz

      fslmaths ${ANATDIR}/gm_tpm_native.nii -thr 0.5 -bin ${ANATDIR}/gm_mask_native.nii
      gunzip -f ${ANATDIR}/gm_mask_native.nii.gz

      3dcalc -a ${ANATDIR}/wm_mask_native.nii -b ${ANATDIR}/csf_mask_native.nii -expr 'a+b' -prefix ${ANATDIR}/nongm_mask_native.nii

      fslmeants -i ${OUTDIR}/proc_data_native${SUFFIX}.nii -m ${ANATDIR}/csf_mask_native.nii -o ${OUTDIR}/csf_sig.txt
      fslmeants -i ${OUTDIR}/proc_data_native${SUFFIX}.nii -m ${ANATDIR}/wm_mask_native.nii -o ${OUTDIR}/wm_sig.txt

      rm -f ${OUTDIR}/nat_mask.nii
      3dAutomask -prefix ${OUTDIR}/nat_mask.nii ${OUTDIR}/proc_data_native.nii
      fslmeants -i ${OUTDIR}/proc_data_native${SUFFIX}.nii -m ${ANATDIR}/nat_mask.nii -o ${OUTDIR}/global_sig.txt

      python run_acompcor.py -d ${OUTDIR} -i ${OUTDIR}/proc_data_native${SUFFIX}.nii -n ${ANATDIR}/nongm_mask_native.nii -b ${OUTDIR}/nat_mask.nii -t ${TR}

      # Remove headers from compcor results
      sed '1d' ${OUTDIR}/acompcor.txt > ${OUTDIR}/tmp_acompcor.txt
      mv ${OUTDIR}/tmp_acompcor.txt ${OUTDIR}/acompcor.txt

      sed '1d' ${OUTDIR}/tcompcor.txt > ${OUTDIR}/tmp_tcompcor.txt
      mv ${OUTDIR}/tmp_tcompcor.txt ${OUTDIR}/tcompcor.txt

      python rp_nuis_calc.py -d ${OUTDIR} -r motion_estimate.par -t ${FD_THR}

    } &> ${LOGDIR}/NuisanceSignal_Calc.log

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
        SUFFIX=''
        if test -f "${OUTDIR}/proc_data_native_nlm.nii"; then
              SUFFIX='_nlm'
              rm -f ${FIXDIR}/filtered_func_data.nii.gz
              cp -f ${OUTDIR}/proc_data_native${SUFFIX}.nii ${FIXDIR}/filtered_func_data.nii
              gzip -f ${FIXDIR}/filtered_func_data.nii
        fi
        CURDIR=`pwd`
        cd ${FIXBIN}
        source ${FIXBIN}/fix -c ${FIXDIR} ${FIX_CLASSIFIER} ${FIX_THR}

        # Non-aggressive clean-up
        source ${FIXBIN}/fix -a ${FIXDIR}/fix4melview_${FIX_CL_LABEL}_thr${FIX_THR}.txt -A
        cd ${CURDIR}

        mv ${FIXDIR}/filtered_func_data_clean.nii.gz ${OUTDIR}/proc_data_native_fix.nii.gz
        gunzip -f ${OUTDIR}/proc_data_native_fix.nii.gz


        python ${QCDIR}/save_image_diff.py -o ${OUTDIR} -i ${OUTDIR} -a proc_data_native.nii -b proc_data_native_fix.nii -msg1 'Before FIX' -msg2 'After FIX'

        #python ${QCDIR}/save_image_diff.py -o ${OUTDIR} -i ${OUTDIR} -a proc_data_native${SUFFIX}.nii -b  proc_data_native_fix.nii -type std -msg1 'Before FIX' -msg2 'After FIX'

        #mv ${OUTDIR}/*.png ${IMGDIR}
        #mv ${OUTDIR}/*.gif ${IMGDIR}
        ## Physiological noise estimation using PHYCAA+ [TODO REF]
        #if [ "${DO_PEST}" -eq "2" ]; then
        #      matlab  "-nodesktop -nosplash " <<<"apply_phycaa_rsn(${TR}, '${OUTDIR}', 'proc_data_native_fix', '${OUTDIR}/nat_mask.nii', '${OUTDIR}/csf_mask_nat.nii' ); exit;"
        #fi
        # Write out the IC time course

        python ${WDIR}/util/save_ics.py -out ${OUTDIR} -class ${FIXDIR}/fix4melview_${FIX_CL_LABEL}_thr${FIX_THR}.txt

      }  &> ${LOGDIR}/FIX_NuisanceCalc.log
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
        python rp_nuis_calc.py -d ${OUTDIR} -r motion_estimate.par -t ${FD_THR}

        rm -f -r ${OUTDIR}/dicer
        rm -f -r ${OUTDIR}/dicer_fix
        CURDIR=`pwd`

        # TODO Set up conf for DiCER PATH
        cd ${DICERPATH}

        SUF=(native_fix native)
        IDXS=0
        mkdir ${OUTDIR}/dicer_fix
        #mkdir ${OUTDIR}/dicer


        INPREF=proc_data_native


        # If background is not 0, fast has problem getting the correct tissue distrub
        #fslmaths ${OUTDIR}/proc_data_native.nii -thr 10  ${OUTDIR}/proc_data_native_tmp.nii
        #gunzip -f ${OUTDIR}/proc_data_native_tmp.nii.gz
        #mv ${OUTDIR}/proc_data_native_tmp.nii ${OUTDIR}/proc_data_native.nii



        cp -f ${OUTDIR}/${INPREF}.nii ${OUTDIR}/tmp.nii
        cp -f ${ANATDIR}/anat_native_brain.nii ${ANATDIR}/anat_native_tmp.nii

        3dresample -prefix ${OUTDIR}/tmp_res.nii -orient lpi -input ${OUTDIR}/tmp.nii
        3dresample -prefix ${ANATDIR}/anat_native_tmp_res.nii -orient lpi -input ${ANATDIR}/anat_native_tmp.nii

        rm -f ${ANATDIR}/anat_native_tmp_res2.nii
        3dUnifize -prefix ${ANATDIR}/anat_native_tmp_res2.nii -input ${ANATDIR}/anat_native_tmp_res.nii

        mv ${ANATDIR}/anat_native_tmp_res2.nii ${ANATDIR}/anat_native_tmp.nii



        mv ${OUTDIR}/tmp_res.nii ${OUTDIR}/tmp.nii


        fslmaths ${OUTDIR}/tmp.nii -Tmean ${OUTDIR}/tmp_mean.nii
        gunzip -f ${OUTDIR}/tmp_mean.nii.gz

        #flirt -in ${ANATDIR}/anat_native_brain.nii -ref ${ANATDIR}/anat_proc.nii  -out ${ANATDIR}/anat_native_brain_fsl \
        #     -usesqform -dof 6 -coarsesearch 18 -finesearch 9 -searchrx -20 20 -searchry -20 20 -searchrz -20 20
        #flirt -applyxfm -init ${ANATDIR}/anat2func_fsl.mat -in ${ANATDIR}/anat_proc_brain.nii -out ${ANATDIR}/anat_native_brain_fsl_func.nii.gz -ref ${ANATDIR}/mean_func_data.nii
        #gunzip -f ${ANATDIR}/anat_native_brain_fsl_func.nii.gz
        #cp -f ${ANATDIR}/anat_native_brain_fsl_func.nii ${ANATDIR}/tmp_anat_native.nii
        #copy_header ${OUTDIR}/tmp_mean.nii ${ANATDIR}/tmp_anat_native.nii
        #rm -f ${OUTDIR}/tmp_mean.nii
        copy_header ${OUTDIR}/tmp_mean.nii ${ANATDIR}/anat_native_tmp.nii
        gzip -f ${OUTDIR}/tmp.nii
        rm -f ${OUTDIR}/tmp.nii

        #cp ${OUTDIR}/non_noise_ics.txt ${OUTDIR}/confounders.txt
        #cp ${OUTDIR}/non_noise_ics.txt ${OUTDIR}/dicer_fix/confounders.txt

        if test -f ${OUTDIR}/dicer_fix/confounders.txt; then
          source DiCER_lightweight.sh -i ${OUTDIR}/tmp.nii.gz -a ${ANATDIR}/anat_native_tmp.nii -w ${OUTDIR}/dicer_fix -s SUBJECT_1_FIX  -p 0 -d -c confounders.txt
        else
          source DiCER_lightweight.sh -i ${OUTDIR}/tmp.nii.gz -a ${ANATDIR}/anat_native_tmp.nii -w ${OUTDIR}/dicer_fix -s SUBJECT_1_FIX  -p 0 -d
        fi

#        gunzip -f /home/luna.kuleuven.be/u0101486/workspace/data/RepImpact/tmp/N1_02/dicer_fix/tmp_detrended_hpf.nii.gz
#        gunzip -f /home/luna.kuleuven.be/u0101486/workspace/data/RepImpact/tmp/N1_02/dicer_fix/tmp_detrended_hpf_dbscan.nii.gz
#
#        python ${QCDIR}/save_image_diff.py -o ${OUTDIR}/dicer_fix -i ${OUTDIR}/dicer_fix -a tmp_detrended_hpf.nii \
#         -b  tmp_detrended_hpf_dbscan.nii \
#         -type std -msg1 'Before DiCER' -msg2 'After DiCER'


        cd ${CURDIR}
        rm ${ANATDIR}/anat_native_tmp.nii


        cp -f ${OUTDIR}/dicer_fix/SUBJECT_1_FIX_dbscan_liberal_regressors.tsv ${OUTDIR}/DiCER_regressors.txt

        rm -f -r ${OUTDIR}/dicer_fix
        rm -f ${OUTDIR}/tmp.nii.gz
        rm -f ${OUTDIR}/tmp.nii
        rm -f ${OUTDIR}/proc_data_native_fix_d.nii

        #python ${QCDIR}/save_image_diff.py -o ${OUTDIR} -i ${OUTDIR} -a proc_data_native_fix.nii -b  proc_data_native_fix_d.nii -type std -msg1 'Before DiCER' -msg2 'After DiCER'
        #mv ${OUTDIR}/*.png ${IMGDIR}
        #mv ${OUTDIR}/*.gif ${IMGDIR}

      }  &> ${LOGDIR}/DiCER.log
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

if [ "$DO_QA" -eq "1" ] || [ "$DO_QA" -eq "2" ]; then
      # Performs QS in native space
      printf "[$SUB] Generating QC plots ... "
      START=$(date -u +%s.%N)


          #fslmaths ${OUTDIR}/proc_data_native_nlm.nii -Tmean ${OUTDIR}/func_mean.nii
          #gunzip -f ${OUTDIR}/func_mean.nii.gz
          #rm -f -r ${OUTDIR}/QA_*
          {
          #cd ${FIXBIN}
          ## Non-aggressive clean-up
          #source ${FIXBIN}/fix -a ${FIXDIR}/fix4melview_${FIX_CL_LABEL}_thr${FIX_THR}.txt -A
          #cd ${CURDIR}

          #rm -f ${OUTDIR}/proc_data_native_fix.nii
          #cp ${FIXDIR}/filtered_func_data_clean.nii.gz ${OUTDIR}/proc_data_native_fix.nii.gz
          #gunzip -f ${OUTDIR}/proc_data_native_fix.nii.gz

          if [ "$DO_QA" -eq "1" ]; then
            rm -f -r ${OUTDIR}/QA_*
          fi
          rm -f ${ANATDIR}/lg400_cobra_native_*.nii

          rm -f ${ANATDIR}/yeo17cobra_native_*.nii

          WDIR=$(awk -F\=  '/^WDIR/{print $2}' "${CONFIG}" | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')


          python ${WDIR}/rp_nuis_calc.py -d ${OUTDIR} -r motion_estimate.par -t 0.3

          model_qa()
          {
              function append_col()
              {
                  FILEA=$1
                  FILEB=$2
                  OUTDIR=$3

                  TMPFILE=tmp_reg_${RANDOM}.tmp


                  pr -mts' ' ${FILEA} ${FILEB} > $TMPFILE

                  mv ${TMPFILE} ${FILEB}

              }

              copy_header()
              {
                  IMBASE=$1
                  IMTARGET=$2

                  # Copy all relevant orientation and position fields
                  QFORM_CODE=`nifti_tool -quiet -disp_hdr -field qform_code -infiles ${IMBASE}`
                  SFORM_CODE=`nifti_tool -quiet -disp_hdr -field sform_code -infiles ${IMBASE}`

                  QUATERNB=`nifti_tool -quiet -disp_hdr -field quatern_b -infiles ${IMBASE}`
                  QUATERNC=`nifti_tool -quiet -disp_hdr -field quatern_c -infiles ${IMBASE}`
                  QUATERND=`nifti_tool -quiet -disp_hdr -field quatern_d -infiles ${IMBASE}`

                  QOFFX=`nifti_tool -quiet -disp_hdr -field qoffset_x -infiles ${IMBASE}`
                  QOFFY=`nifti_tool -quiet -disp_hdr -field qoffset_y -infiles ${IMBASE}`
                  QOFFZ=`nifti_tool -quiet -disp_hdr -field qoffset_z -infiles ${IMBASE}`

                  SROWX=`nifti_tool -quiet -disp_hdr -field srow_x -infiles ${IMBASE}`
                  SROWY=`nifti_tool -quiet -disp_hdr -field srow_y -infiles ${IMBASE}`
                  SROWZ=`nifti_tool -quiet -disp_hdr -field srow_z -infiles ${IMBASE}`

                  PIXDIM=`nifti_tool -quiet -disp_hdr -field pixdim -infiles ${IMBASE}`

                  SLSTART=`nifti_tool -quiet -disp_hdr -field slice_start -infiles ${IMBASE}`
                  SLEND=`nifti_tool -quiet -disp_hdr -field slice_end -infiles ${IMBASE}`
                  SLCODE=`nifti_tool -quiet -disp_hdr -field slice_code -infiles ${IMBASE}`

                  nifti_tool -mod_hdr -overwrite -infiles ${IMTARGET}  \
                        -mod_field qform_code ${QFORM_CODE} -mod_field sform_code ${SFORM_CODE} \
                        -mod_field quatern_b ${QUATERNB} -mod_field quatern_c ${QUATERNC} \
                        -mod_field quatern_d ${QUATERND} -mod_field qoffset_x ${QOFFX} \
                        -mod_field qoffset_y ${QOFFY} -mod_field qoffset_z ${QOFFZ} \
                        -mod_field srow_x "${SROWX}" -mod_field srow_y "${SROWY}" \
                        -mod_field srow_z "${SROWZ}" -mod_field pixdim "${PIXDIM}" \
                        -mod_field slice_start ${SLSTART} -mod_field slice_end ${SLEND} \
                        -mod_field slice_code ${SLCODE}
              }

              OUTDIR=${1}
              ANATDIR=${2}
              DARTEL_TEMPLATE_PREF=${3}
              WDIR=${4}
              TR=${5}
              QCDIR=${6}
              INPREF=${7}
              REG_AGG=${8}
              REG_MODEL=${9}
              DELMOD=${10}

              QA_DIR=${OUTDIR}/QA_${REG_MODEL}_${REG_AGG}
              mkdir -p ${QA_DIR}/


              rm -f ${OUTDIR}/${INPREF}_${REG_MODEL}.nii
              rm -f ${OUTDIR}/${INPREF}_${REG_MODEL}.nii.gz
              if [ "$DELMOD" -eq "1" ]; then
                rm -f ${QA_DIR}/*
              fi


              rm -f ${OUTDIR}/regressors_${REG_MODEL}.txt
              touch ${OUTDIR}/regressors_${REG_MODEL}.txt


              if [[ ${REG_MODEL} == *"RP6"* ]]; then
                  append_col ${OUTDIR}/motion_estimate.par ${OUTDIR}/regressors_${REG_MODEL}.txt ${OUTDIR}
              fi

              if [[ ${REG_MODEL} == *"RP24"* ]]; then
                  append_col ${OUTDIR}/motion_regressors_12.txt  ${OUTDIR}/regressors_${REG_MODEL}.txt ${OUTDIR}
                  append_col ${OUTDIR}/motion_regressors_24.txt  ${OUTDIR}/regressors_${REG_MODEL}.txt ${OUTDIR}
              fi

              if [[ ${REG_MODEL} == *"CC"* ]] || [[ ${REG_MODEL} == *"DD"* ]]; then
                  3dcalc -a ${ANATDIR}/wm_tpm_native.nii -b ${ANATDIR}/csf_tpm_native.nii -expr 'astep(a+b, 0.95)' -prefix ${ANATDIR}/nongm_mask_native_${REG_MODEL}.nii
                  #cp ${ANATDIR}/nongm_mask_native.nii ${ANATDIR}/nongm_mask_native_${REG_MODEL}.nii
                  copy_header ${OUTDIR}/${INPREF} ${ANATDIR}/nongm_mask_native_${REG_MODEL}.nii

                  if [[ ${REG_MODEL} == *"CC"* ]]; then
                        python ${WDIR}/run_acompcor.py -d ${OUTDIR} -i ${OUTDIR}/${INPREF}.nii -n ${ANATDIR}/nongm_mask_native_${REG_MODEL}.nii -b ${OUTDIR}/nat_mask.nii -t ${TR} \
                                      -var 0.25 -aout acompcor_${REG_MODEL} -tout tcompcor_${REG_MODEL}
                  else
                        python ${WDIR}/run_acompcor.py -d ${OUTDIR} -i ${OUTDIR}/${INPREF}.nii -n ${ANATDIR}/nongm_mask_native_${REG_MODEL}.nii -b ${OUTDIR}/nat_mask.nii -t ${TR} \
                                      -var 0.5 -aout acompcor_${REG_MODEL} -tout tcompcor_${REG_MODEL}

                  fi

                  sed '1d' ${OUTDIR}/acompcor_${REG_MODEL}.txt > ${OUTDIR}/tmp_acompcor_${REG_MODEL}.txt
                  mv ${OUTDIR}/tmp_acompcor_${REG_MODEL}.txt ${OUTDIR}/acompcor_${REG_MODEL}.txt

                  sed '1d' ${OUTDIR}/tcompcor_${REG_MODEL}.txt > ${OUTDIR}/tmp_tcompcor_${REG_MODEL}.txt
                  mv ${OUTDIR}/tmp_tcompcor_${REG_MODEL}.txt ${OUTDIR}/tcompcor_${REG_MODEL}.txt

                  #append_col ${OUTDIR}/acompcor.txt  ${OUTDIR}/regressors_${REG_MODEL}.txt ${OUTDIR}
                  #append_col ${OUTDIR}/tcompcor_${REG_MODEL}.txt  ${OUTDIR}/regressors_${REG_MODEL}.txt ${OUTDIR}
                  if [[ ${REG_MODEL} == *"CC"* ]]; then
                    append_col ${OUTDIR}/acompcor_${REG_MODEL}.txt  ${OUTDIR}/regressors_${REG_MODEL}.txt ${OUTDIR}
                  else
                    append_col ${OUTDIR}/tcompcor_${REG_MODEL}.txt  ${OUTDIR}/regressors_${REG_MODEL}.txt ${OUTDIR}
                  fi
                  rm -f ${ANATDIR}/nongm_mask_native_${REG_MODEL}.nii
              fi

              if [[ ${REG_MODEL} == *"WM"* ]]; then
                  fslmaths ${ANATDIR}/wm_tpm_native.nii -thr 0.9 -bin ${ANATDIR}/wm_mask_native.nii
                  gunzip -f ${ANATDIR}/wm_mask_native.nii.gz

                  fslmeants -i ${OUTDIR}/${INPREF}.nii -m ${ANATDIR}/wm_mask_native.nii -o ${OUTDIR}/wm_sig_${REG_MODEL}.txt

                  append_col ${OUTDIR}/wm_sig_${REG_MODEL}.txt  ${OUTDIR}/regressors_${REG_MODEL}.txt ${OUTDIR}
                  rm -f ${OUTDIR}/wm_sig_${REG_MODEL}.txt
              fi


              if [[ ${REG_MODEL} == *"CSF"* ]]; then
                  fslmaths ${ANATDIR}/csf_tpm_native.nii -thr 0.9 -bin ${ANATDIR}/csf_mask_native.nii
                  gunzip -f ${ANATDIR}/csf_mask_native.nii.gz

                  fslmeants -i ${OUTDIR}/${INPREF}.nii -m ${ANATDIR}/csf_mask_native.nii -o ${OUTDIR}/csf_sig_${REG_MODEL}.txt

                  append_col ${OUTDIR}/csf_sig_${REG_MODEL}.txt  ${OUTDIR}/regressors_${REG_MODEL}.txt ${OUTDIR}
                  rm -f ${OUTDIR}/csf_sig_${REG_MODEL}.txt
              fi

              #if [[ ${REG_MODEL} == *"FIX"* ]]; then
              #    append_col ${OUTDIR}/noise_ics.txt  ${OUTDIR}/regressors_${REG_MODEL}.txt ${OUTDIR}
              #fi

              if [[ ${REG_MODEL} == *"DiC"* ]]; then
                  append_col ${OUTDIR}/DiCER_regressors.txt  ${OUTDIR}/regressors_${REG_MODEL}.txt ${OUTDIR}
              fi

              if [[ ${REG_MODEL} == *"FAC"* ]]; then
                  append_col ${OUTDIR}/noise_ics.txt  ${OUTDIR}/regressors_${REG_MODEL}.txt ${OUTDIR}
              fi

              if [[ ${REG_MODEL} == *"FAD"* ]]; then

                python ${WDIR}/util/orthogonalize_regressors.py \
                        -regressors ${OUTDIR}/noise_ics.txt \
                        -signal ${OUTDIR}/non_noise_ics.txt \
                        -out ${OUTDIR}/noise_ics_${REG_MODEL}.txt

                  append_col ${OUTDIR}/noise_ics_${REG_MODEL}.txt  ${OUTDIR}/regressors_${REG_MODEL}.txt ${OUTDIR}
                  rm -f ${OUTDIR}/noise_ics_${REG_MODEL}.txt
              fi

              if [[ ${REG_MODEL} == *"GSR"* ]]; then
                  append_col ${OUTDIR}/global_sig.txt  ${OUTDIR}/regressors_${REG_MODEL}.txt ${OUTDIR}
              fi

              if [[ ${REG_MODEL} == *"PHY"* ]]; then
                  append_col ${OUTDIR}/physio.txt  ${OUTDIR}/regressors_${REG_MODEL}.txt ${OUTDIR}
              fi

              3dcalc -a ${ANATDIR}/gm_tpm_native.nii  -expr 'astep(a, 0.3)' -prefix ${ANATDIR}/gm_mask_native_${REG_MODEL}.nii
              copy_header ${OUTDIR}/${INPREF} ${ANATDIR}/gm_mask_native_${REG_MODEL}.nii


              #if [[ ${REG_MODEL} == "NONE" ]]; then
              CMD="3dTproject -input ${OUTDIR}/${INPREF}.nii -TR ${TR} "
              #else
              #    CMD="3dTproject -input ${OUTDIR}/${INPREF}_${REG_MODEL}.nii -TR ${TR} -norm "
              #fi

              CMD="${CMD} -polort 2 "

              if [ "$DELMOD" -eq "1" ]; then
                  CMD="${CMD} -mask ${ANATDIR}/gm_mask_native_${REG_MODEL}.nii "
              else
                  CMD="${CMD} -mask ${OUTDIR}/nat_mask_tpm.nii "
              fi
              CMD="${CMD} -prefix ${OUTDIR}/${INPREF}_${REG_MODEL}.nii "




              if [[ ${REG_MODEL} == *"CEN"* ]]; then
                  CMD="${CMD} -cenmode KILL -censor ${OUTDIR}/temporal_mask_fd.txt "
              fi

              if [[ ${REG_MODEL} == *"PCA"* ]]; then
                 #3dpc -prefix ${OUTDIR}/regressors_${REG_MODEL}_PCA -eigonly -pcsave ALL -dmean -nscale ${OUTDIR}/regressors_${REG_MODEL}.txt
                 if [[ ${REG_MODEL} == *"PCA80"* ]]; then
                   python ${WDIR}/util/nusiance_pca.py -outfile ${OUTDIR}/regressors_${REG_MODEL}_PCA.txt -varexp 80 -nuis_file ${OUTDIR}/regressors_${REG_MODEL}.txt
                 fi

                 if [[ ${REG_MODEL} == *"PCA95"* ]]; then
                   python ${WDIR}/util/nusiance_pca.py -outfile ${OUTDIR}/regressors_${REG_MODEL}_PCA.txt -varexp 95 -nuis_file ${OUTDIR}/regressors_${REG_MODEL}.txt
                 fi

                 #rm -f ${OUTDIR}/regressors_${REG_MODEL}_PCA_eig.1D
                 cp ${OUTDIR}/regressors_${REG_MODEL}_PCA.txt ${OUTDIR}/regressors_${REG_MODEL}.txt
              fi

              #CMD="${CMD} -ort ${OUTDIR}/regressors_${REG_MODEL}_ortho.txt "


              if [[ ${REG_MODEL} == *"DVOX"* ]]; then
                    CMD="${CMD} -dsort ${OUTDIR}/voxel_displ_DX.nii "
                    CMD="${CMD} -dsort ${OUTDIR}/voxel_displ_DY.nii "
                    CMD="${CMD} -dsort ${OUTDIR}/voxel_displ_DZ.nii "
              fi

              echo "${REG_MODEL}"
              echo "RUNNING: ${CMD}"
              if [ -s "${OUTDIR}/regressors_${REG_MODEL}.txt" ]; then
                  CMD="${CMD} -ort ${OUTDIR}/regressors_${REG_MODEL}.txt "
              fi


              eval "${CMD}"
              rm -f ${OUTDIR}/regressors_${REG_MODEL}.txt
              rm -f ${OUTDIR}/tcompcor_${REG_MODEL}.txt
              rm -f ${OUTDIR}/acompcor_${REG_MODEL}.txt
              rm -f ${OUTDIR}/tcompcor_${REG_MODEL}_meta.txt
              rm -f ${OUTDIR}/acompcor_${REG_MODEL}_meta.txt
              3dDespike -nomask -NEW -localedit -cut 1.96 2.96 \
                    -prefix  ${OUTDIR}/${INPREF}_${REG_MODEL}_tmp.nii \
                            ${OUTDIR}/${INPREF}_${REG_MODEL}.nii

              mv ${OUTDIR}/${INPREF}_${REG_MODEL}_tmp.nii ${OUTDIR}/${INPREF}_${REG_MODEL}.nii

              copy_header ${OUTDIR}/${INPREF} ${OUTDIR}/${INPREF}_${REG_MODEL}.nii

              3dTto1D -input ${OUTDIR}/${INPREF}_${REG_MODEL}.nii -method DVARS \
                  -mask ${ANATDIR}/gm_mask_native_${REG_MODEL}.nii \
                  -prefix ${QA_DIR}/dvar.1d

              rm -f ${ANATDIR}/gm_mask_native_${REG_MODEL}.nii
              #rm -f ${OUTDIR}/csf_sig_${REG_MODEL}.txt
              #rm -f ${OUTDIR}/wm_sig_${REG_MODEL}.txt


              #NVOLS=`fslnvols ${OUTDIR}/${INPREF}_${REG_MODEL}.nii`
              #matlab "-nodesktop -nosplash " <<<"coreg_normalise('${ANATDIR}/mean_func_norm.nii', '${OUTDIR}/${INPREF}_${REG_MODEL}.nii', ${NVOLS}, '${ANATDIR}/anat_proc.nii', [0 0 1], {'${DARTEL_TEMPLATE_PREF}1.nii', '${DARTEL_TEMPLATE_PREF}2.nii', '${DARTEL_TEMPLATE_PREF}3.nii', '${DARTEL_TEMPLATE_PREF}4.nii', '${DARTEL_TEMPLATE_PREF}5.nii', '${DARTEL_TEMPLATE_PREF}6.nii'}, [5 5 5]); exit;"
              if [ "$DELMOD" -eq "2" ]; then
                3dReHo  -inset ${OUTDIR}/${INPREF}_${REG_MODEL}.nii  -prefix ${QA_DIR}/REHO.nii
                copy_header ${ANATDIR}/mean_func_data_nds.nii ${QA_DIR}/REHO.nii
                matlab "-nodesktop -nosplash " <<<"coreg_normalise('${ANATDIR}/mean_func_data_nds.nii', '${QA_DIR}/REHO.nii', 1, '${ANATDIR}/anat_proc.nii', [0 0 1], {'${DARTEL_TEMPLATE_PREF}1.nii', '${DARTEL_TEMPLATE_PREF}2.nii', '${DARTEL_TEMPLATE_PREF}3.nii', '${DARTEL_TEMPLATE_PREF}4.nii', '${DARTEL_TEMPLATE_PREF}5.nii', '${DARTEL_TEMPLATE_PREF}6.nii'}, [6 6 6]); exit;"
                mv ${QA_DIR}/swREHO.nii ${QA_DIR}/REHO.nii
              fi

              #rm -f ${ANATDIR}/w${INPREF}_${REG_MODEL}.nii

            rm -f  ${ANATDIR}/func_${REG_MODEL}.nii
            rm -f ${ANATDIR}/lg400_cobra_native_${REG_MODEL}.nii

            cp ${ANATDIR}/lg400_cobra_native.nii ${ANATDIR}/lg400_cobra_native_${REG_MODEL}.nii
            cp ${ANATDIR}/yeo17cobra_native.nii ${ANATDIR}/yeo17cobra_native_${REG_MODEL}.nii
            #3dTcat -prefix ${ANATDIR}/func_${REG_MODEL}.nii  ${OUTDIR}/${INPREF}_${REG_MODEL}.nii[0]

            copy_header ${OUTDIR}/${INPREF}_${REG_MODEL}.nii ${ANATDIR}/lg400_cobra_native_${REG_MODEL}.nii
            copy_header ${OUTDIR}/${INPREF}_${REG_MODEL}.nii ${ANATDIR}/yeo17cobra_native_${REG_MODEL}.nii


            #fslmaths ${OUTDIR}/${INPREF}_${REG_MODEL}.nii -nan ${OUTDIR}/${INPREF}_${REG_MODEL}_b.nii
            #gunzip -f ${OUTDIR}/${INPREF}_${REG_MODEL}_b.nii.gz
            #mv ${OUTDIR}/${INPREF}_${REG_MODEL}_b.nii ${OUTDIR}/${INPREF}_${REG_MODEL}.nii
            #python ${QCDIR}/save_image_diff.py -o ${QA_DIR} -i ${OUTDIR} -a ${INPREF}.nii -b ${INPREF}_${REG_MODEL}.nii -msg1 "Raw" -msg2 "${REG_MODEL}"

            if [ "$DELMOD" -eq "2" ]; then
                  python ${QCDIR}/QC_grey_plot.py -out ${QA_DIR} -in ${OUTDIR}  -mpe motion_estimate.par \
                    -fname ${INPREF}_${REG_MODEL}.nii -csf_name ${ANATDIR}/csf_mask_group.nii \
                    -wm_name ${ANATDIR}/wm_mask_group.nii \
                    -gm_name ${ANATDIR}/gm_mask_group.nii \
                    -norm zscore -range 1.0 \
                    -prog AFNI -outf 02_Greyplot_${REG_MODEL}  -dpi 100 -tr ${TR}
            fi

            if [ "$DELMOD" -eq "1" ]; then
              python ${QCDIR}/QC_FC.py -out ${QA_DIR} -in ${OUTDIR}  \
                          -fname ${INPREF}_${REG_MODEL}.nii -outf FC_${REG_MODEL}  \
                          -atlas ${ANATDIR}/lg400_cobra_native_${REG_MODEL}.nii \
                          -atlas_name local_global_cobra -TR ${TR} \
                          -dist /home/luna.kuleuven.be/u0101486/workspace/fmri_proc/atlas/CAT12/lg400_cobra_distance.txt \
                          -vmin -0.6 -vmax 0.6 \
                           -dpi 72

            fi



     #        python ${QCDIR}/QC_FC.py -out ${QA_DIR} -in ${OUTDIR}  \
    #                    -fname ${INPREF}_${REG_MODEL}.nii -outf FC_${REG_MODEL}  \
   #                     -atlas ${ANATDIR}/yeo17cobra_native_${REG_MODEL}.nii \
  #                      -atlas_name yeo17_cobra \
 #                       -vmin -0.6 -vmax 0.6 \
#                         -dpi 72

              if [ "$DELMOD" -eq "1" ]; then
                rm -f ${OUTDIR}/${INPREF}_${REG_MODEL}.nii
              fi
              rm -f  ${ANATDIR}/func_${REG_MODEL}.nii

              rm -f ${ANATDIR}/lg400_cobra_native_${REG_MODEL}.nii
              rm -f ${ANATDIR}/yeo17cobra_native_${REG_MODEL}.nii
              rm -f ${OUTDIR}/${INPREF}_${REG_MODEL}.nii.gz

              rm -f ${OUTDIR}/regressors_${REG_MODEL}_ortho.txt

              rm -f ${OUTDIR}/sw${INPREF}_${REG_MODEL}.nii
              rm -f ${OUTDIR}/sw${INPREF}_${REG_MODEL}.nii.gz
      }
      export -f model_qa
#NONAGG
      #parallel -j4 --line-buffer model_qa ::: ${OUTDIR} ::: ${ANATDIR} ::: ${DARTEL_TEMPLATE_PREF} ::: QA_FIX ::: ${TR} ::: ${QCDIR}  :::  proc_data_mni_fix ::: SFIX SFIX_D SRP24WM1CSF1 SRP24WM1CSF1_D SRP24CC SRP9
      #parallel -j4 --line-buffer model_qa ::: ${OUTDIR} ::: ${ANATDIR} ::: ${DARTEL_TEMPLATE_PREF} ::: ${WDIR} ::: ${TR} ::: ${QCDIR}  :::  proc_data_native_fix ::: AGG  ::: FIX FIX_CC FIX_DiC FIX_RP24 FIX_DiC_RP24 FIX_CC_RP24
      #FAD_DiC FAD_DiC_DVOX FAD_DiC_DVOX_RP6 FAD_DiC_RP24 FAD_PHY FAD_PHY_DVOX FAD_PHY_DVOX_RP6 FAD_PHY_RP24
      if [ "$DO_QA" -eq "1" ]; then
          if test -f ${OUTDIR}/proc_data_native_u.nii; then
              parallel -j4 --line-buffer model_qa ::: ${OUTDIR} ::: ${ANATDIR} ::: ${DARTEL_TEMPLATE_PREF} ::: ${WDIR} ::: ${TR} ::: ${QCDIR}  :::  proc_data_native_u :::   AGG ::: FAC_WM_CSF FAC_DD FAD_DiC_RP24 FAC_DiC_RP24 FAD_DD_RP24 FAC_DD_RP24 ::: ${DO_QA}
          else
              # FAC_CC_RP6 FAD_CC_RP6 FAC_DiC_RP6 FAD_DiC_RP6 FAC_CC_RP24 FAD_CC_RP24 FAC_DiC_RP24 FAD_DiC_RP24 FAC_WM_CSF_RP6 FAD_WM_CSF_RP24 RP24_WM_CSF RP24_CC FAC_PHY FAD_PHY RP24_PHY
              parallel -j1 --line-buffer model_qa ::: ${OUTDIR} ::: ${ANATDIR} ::: ${DARTEL_TEMPLATE_PREF} ::: ${WDIR} ::: ${TR} ::: ${QCDIR}  :::  proc_data_native ::: AGG ::: \
                  RP24_WM_CSF_CEN FAD_DiC_RP24_CEN FAC_DiC_RP24_CEN \
                  FAD_DiC_RP6_CEN FAC_DiC_RP6_CEN \
                  FAD_DiC_RP6_DVOX_CEN FAC_DiC_RP6_DVOX_CEN \
                  FAD_DiC_RP24_CEN_DVOX FAC_DiC_RP24_CEN_DVOX \
                  FAD_CC_RP24_CEN_DVOX FAC_CC_RP24_CEN_DVOX \
                  FAD_DiC_CC_RP24_PCA95_CEN_DVOX  \
                  FAC_DiC_CC_RP24_PCA95_CEN_DVOX \
                  ::: ${DO_QA}
          fi
      else
          parallel -j1 --line-buffer model_qa ::: ${OUTDIR} ::: ${ANATDIR} ::: ${DARTEL_TEMPLATE_PREF} ::: ${WDIR} ::: ${TR} ::: ${QCDIR}  :::  proc_data_native_u :::  AGG ::: ${REG_MODEL} ::: ${DO_QA}
      fi
      rm -f ${OUTDIR}/proc_data_native_fix.nii
      rm -f ${OUTDIR}/ptmp_*.nii
      rm -f ${OUTDIR}/tmp_*.nii
} &> ${OUTDIR}/QA.log


      END=$(date -u +%s.%N)
      DIFF=`echo "( $END - $START )" | bc`
      printf "DONE [%.1f s]\n" $DIFF
fi



# Apply normalisation WARP
if [ "$DO_NORM" -eq "1" ]; then

      printf "[$SUB] Warping functional images to MNI space ... "
      START=$(date -u +%s.%N)

      {
            WDIR=$(awk -F\=  '/^WDIR/{print $2}' "${CONFIG}" | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
            cd ${WDIR}

            rm -f ${ANATDIR}/c6anat_proc.nii

            rm -f ${OUTDIR}/proc_data_mni*.nii


            INPREF=proc_data_native_u_${REG_MODEL}

            rm -f ${OUTDIR}/fix_mean.nii

            cp -f ${OUTDIR}/${INPREF}.nii ${OUTDIR}/proc_data_native_fix_tmp.nii
            NVOLS=`fslnvols ${OUTDIR}/proc_data_native_fix_tmp.nii`
            3dTcat -prefix ${OUTDIR}/fix_mean.nii ${OUTDIR}/proc_data_native_fix_tmp.nii[0]
            cp -f ${ANATDIR}/mean_func_data_nds.nii ${OUTDIR}/mean_func_data_nds.nii

            #matlab "-nodesktop -nosplash " <<<"coreg_same_image('${OUTDIR}/mean_func_data_nds.nii', '${OUTDIR}/fix_mean.nii', '${OUTDIR}/proc_data_native_fix_tmp.nii', ${NVOLS}); exit;"
            matlab "-nodesktop -nosplash " <<<"coreg_normalise('${ANATDIR}/mean_func_data_nds.nii', {'${OUTDIR}/proc_data_native_fix_tmp.nii'}, ${NVOLS}, '${ANATDIR}/anat_proc.nii', [0 0 1], {'${DARTEL_TEMPLATE_PREF}1.nii', '${DARTEL_TEMPLATE_PREF}2.nii', '${DARTEL_TEMPLATE_PREF}3.nii', '${DARTEL_TEMPLATE_PREF}4.nii', '${DARTEL_TEMPLATE_PREF}5.nii', '${DARTEL_TEMPLATE_PREF}6.nii'}, [5 5 5], 4); exit;"

            matlab "-nodesktop -nosplash " <<<"coreg_normalise('${ANATDIR}/mean_func_data_nds.nii', {'${ANATDIR}/csf_tpm_native.nii', '${ANATDIR}/wm_tpm_native.nii', '${ANATDIR}/gm_tpm_native.nii'}, ${NVOLS}, '${ANATDIR}/anat_proc.nii', [0 0 1], {'${DARTEL_TEMPLATE_PREF}1.nii', '${DARTEL_TEMPLATE_PREF}2.nii', '${DARTEL_TEMPLATE_PREF}3.nii', '${DARTEL_TEMPLATE_PREF}4.nii', '${DARTEL_TEMPLATE_PREF}5.nii', '${DARTEL_TEMPLATE_PREF}6.nii'}, [2 2 2],1); exit;"
            #matlab "-nodesktop -nosplash " <<<"coreg_normalise('${OUTDIR}/mean_func_data_nds.nii', '${ANATDIR}/wm_tpm_native.nii', ${NVOLS}, '${ANATDIR}/anat_proc.nii', [0 0 1], {'${DARTEL_TEMPLATE_PREF}1.nii', '${DARTEL_TEMPLATE_PREF}2.nii', '${DARTEL_TEMPLATE_PREF}3.nii', '${DARTEL_TEMPLATE_PREF}4.nii', '${DARTEL_TEMPLATE_PREF}5.nii', '${DARTEL_TEMPLATE_PREF}6.nii'}, [0 0 0]); exit;"
            #matlab "-nodesktop -nosplash " <<<"coreg_normalise('${OUTDIR}/mean_func_data_nds.nii', '${ANATDIR}/gm_tpm_native.nii', ${NVOLS}, '${ANATDIR}/anat_proc.nii', [0 0 1], {'${DARTEL_TEMPLATE_PREF}1.nii', '${DARTEL_TEMPLATE_PREF}2.nii', '${DARTEL_TEMPLATE_PREF}3.nii', '${DARTEL_TEMPLATE_PREF}4.nii', '${DARTEL_TEMPLATE_PREF}5.nii', '${DARTEL_TEMPLATE_PREF}6.nii'}, [0 0 0]); exit;"

            rm -f ${OUTDIR}/fix_mean.nii
            rm -f ${OUTDIR}/rfix_mean.nii
            rm -f ${OUTDIR}/proc_data_native_fix_tmp.nii
            rm -f ${OUTDIR}/rproc_data_native_fix_tmp.nii
            rm -f ${OUTDIR}/rproc_data_native_fix_tmp.mat
            rm -f ${OUTDIR}/rcsf_tpm_native.nii
            rm -f ${OUTDIR}/rcsf_tpm_native.mat
            rm -f ${OUTDIR}/rwm_tpm_native.nii
            rm -f ${OUTDIR}/rwm_tpm_native.mat
            rm -f ${OUTDIR}/rgm_tpm_native.nii
            rm -f ${OUTDIR}/rgm_tpm_native.mat

            rm -f ${OUTDIR}/proc_data_native_fix_tmp.mat
            rm -f ${OUTDIR}/mean_func_data_nds.nii
            mv ${OUTDIR}/swproc_data_native_fix_tmp.nii ${OUTDIR}/proc_data_mni.nii
            mv ${ANATDIR}/swcsf_tpm_native.nii ${ANATDIR}/csf_tpm_group.nii
            mv ${ANATDIR}/swwm_tpm_native.nii ${ANATDIR}/wm_tpm_group.nii
            mv ${ANATDIR}/swgm_tpm_native.nii ${ANATDIR}/gm_tpm_group.nii


            #fslmaths ${ANATDIR}/csf_mask_group.nii -thr 0.3 -bin -ero  ${ANATDIR}/csf_mask_group_.nii
            #gunzip -f ${ANATDIR}/csf_mask_group_.nii.gz
            #mv ${ANATDIR}/csf_mask_group_.nii  ${ANATDIR}/csf_mask_group.nii

            #fslmaths ${ANATDIR}/wm_mask_group.nii -thr 0.4 -bin  -ero  ${ANATDIR}/wm_mask_group_.nii
            #gunzip -f ${ANATDIR}/wm_mask_group_.nii.gz
            #mv ${ANATDIR}/wm_mask_group_.nii  ${ANATDIR}/wm_mask_group.nii

            #fslmaths ${ANATDIR}/gm_mask_group.nii -thr 0.5 -bin   ${ANATDIR}/gm_mask_group_.nii
            #gunzip -f ${ANATDIR}/gm_mask_group_.nii.gz
            #mv ${ANATDIR}/gm_mask_group_.nii  ${ANATDIR}/gm_mask_group.nii

            3dAutomask -prefix ${OUTDIR}/mask_mni.nii ${OUTDIR}/proc_data_mni.nii
            #3dBlurToFWHM -quiet -prefix ${OUTDIR}/proc_data_mni_s.nii -input ${OUTDIR}/proc_data_mni.nii -FWHM 6 -rate 2 -mask ${OUTDIR}/mask_mni.nii

            #mv ${OUTDIR}/proc_data_mni_s.nii ${OUTDIR}/proc_data_mni.nii
            rm -f ${OUTDIR}/${INPREF}.nii
            #fslmaths  ${OUTDIR}/proc_data_mni.nii -ing 1000 ${OUTDIR}/proc_data_mni.nii -odt short
            #gunzip -f ${OUTDIR}/proc_data_mni.nii.gz
            #mv ${OUTDIR}/proc_data_mni.nii ${OUTDIR}/proc_data_${CURR_SUFF}.nii



} &> ${OUTDIR}/05_norm2mni.log

      END=$(date -u +%s.%N)
      DIFF=`echo "( $END - $START )" | bc`
      printf "DONE [%.1f s]\n" $DIFF
fi


exit



#rm -f ${OUTDIR}/proc_data_native*.nii
exit

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
