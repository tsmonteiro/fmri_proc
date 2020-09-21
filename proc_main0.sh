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


ANATDIR=${OUTDIR}/anat
IMGDIR=${OUTDIR}/figures
FIXDIR=${OUTDIR}/FIX

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

      mkdir ${IMGDIR}
      mkdir ${ANATDIR}

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


                  3dresample -prefix ${OUTDIR}/tmp.nii -orient ras -input ${OUTDIR}/${PREF}func_data.nii
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

            cp  -f ${OUTDIR}/${PREF}func_data.nii  ${OUTDIR}/m${PREF}func_data.nii
            #3dTcat -prefix ${OUTDIR}/m${PREF}func_data.nii -tr ${TR} ${OUTDIR}/tmp/func_data_n*.nii
            PREF=m${PREF}
            #3dTcat -prefix ${OUTDIR}/${PREF}func_data.nii -tr ${TR} ${OUTDIR}/tmp/func_data_n*.nii


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
                  3dmerge -prefix ${OUTDIR}/ref_vol_smooth.nii -1blur_fwhm 4 ${OUTDIR}/ref_vol.nii


                  3dvolreg -heptic -prefix ${OUTDIR}/r${PREF}func_data.nii -base ${OUTDIR}/ref_vol_smooth.nii -rot_thresh 0.001 -delta 0.15 \
                        -x_thresh 0.001 -zpad 10 -maxite 75 -1Dfile ${OUTDIR}/motion_estimate.par -maxdisp1D ${OUTDIR}/maximum_disp.1d  \
                         ${OUTDIR}/${PREF}func_data.nii



                  1dplot -thick -volreg -png ${OUTDIR}/motion_estimate_0.png -one ${OUTDIR}/motion_estimate.par
                  1dplot -thick -png ${OUTDIR}/maximum_disp_0.png -one ${OUTDIR}/maximum_disp.1d
                  1dplot -thick -png ${OUTDIR}/maximum_disp_delt_0.png -one ${OUTDIR}/maximum_disp.1d_delt

                  mv ${OUTDIR}/*.png ${IMGDIR}
                  mv ${OUTDIR}/*.gif ${IMGDIR}

                  PREF=r${PREF}

                  log_command_div 'MOTION CORRECTION'

            fi


            if [ "${DO_PEST}" == "3" ]; then
                  # At the moment, only working for one specific data set
                  # 4 --> Number of dummy scans. Pass this to parameters
                  matlab "-nodesktop -nosplash " <<<"cd '${WDIR}/matlab'; process_physio_data('${OUTDIR}', $NVOLS, 4, ${TR}, ${NZ}); exit;"

                  mv ${OUTDIR}/slice_timing_g.txt ${OUTDIR}/slice_acq.txt
                  python ./util/write_slice_timing.py -tr ${TR} -nsl ${NSLICES} -acq 3 -out ${OUTDIR}/slice_acq.txt


                  3dREMLfit -input ${OUTDIR}/${PREF}func_data.nii \
                      -addbase ${OUTDIR}/physio_reg.ortho \
                      -GOFORIT \
                      -Rwherr ${OUTDIR}/p${PREF}func_data.nii

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

                  #FIXME slice timing name is now not consistent
                  print_debug 'Slice timing correction' "filtershift --in=${OUTDIR}/${PREF}func_data.nii \
                    --cf=${CF} --rt=${REF_TIME} --TR=${TR} --timing=${OUTDIR}/slice_timing_gs.txt --out=${OUTDIR}/a${PREF}func_data.nii"

                  #--rt=${REF_TIME}
                  filtershift --in=${OUTDIR}/${PREF}func_data.nii --timing=${OUTDIR}/slice_timing_gs.txt --cf=${CF} --rs=30 --TR=${TR} --out=${OUTDIR}/a${PREF}func_data.nii
                  gunzip -f ${OUTDIR}/a${PREF}func_data.nii.gz
                  python ${QCDIR}/save_image_diff.py -o ${OUTDIR} -i ${OUTDIR} -a ${PREF}func_data.nii -b a${PREF}func_data.nii -msg1 'Before STC' -msg2 'After STC'
                  PREF=a${PREF}

                  mv ${OUTDIR}/*.png ${IMGDIR}
                  mv ${OUTDIR}/*.gif ${IMGDIR}
                  log_command_div 'Slice Timing correction'
            fi

            if [ "$DO_DPK" -eq "2" ]; then
                  print_debug 'Despiking [-localedit -NEW]'
                  3dDespike -nomask -NEW -localedit -cut 2 4 -prefix ${OUTDIR}/d${PREF}func_data.nii ${OUTDIR}/${PREF}func_data.nii
                  python ${QCDIR}/save_image_diff.py -o ${OUTDIR} -i ${OUTDIR} -a ${PREF}func_data.nii -b d${PREF}func_data.nii -msg1 'Before Despike' -msg2 'After Despike'
                  PREF=d${PREF}
                  mv ${OUTDIR}/*.png ${IMGDIR}
                  mv ${OUTDIR}/*.gif ${IMGDIR}
                  log_command_div 'Despike'
            fi

            if [ "$DO_DEOBL" -eq "1" ]; then
                  print_debug 'Deobliquing volumes'

                  #TODO CHECK AND CORRECT ORIENTATION INFO!

                  #3dWarp -deoblique -newgrid ${DEOBL_VOX} -NN -prefix ${OUTDIR}/w${PREF}func_data.nii ${OUTDIR}/${PREF}func_data.nii
                  #cp ${OUTDIR}/${PREF}func_data.nii ${OUTDIR}/w${PREF}func_data.nii
                  3dresample -prefix ${OUTDIR}/w${PREF}func_data.nii -orient ras -input ${OUTDIR}/${PREF}func_data.nii
                  3drefit -deoblique ${OUTDIR}/w${PREF}func_data.nii

                  PREF=w${PREF}
                  log_command_div 'DEOBLIQUE'
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

                  3dresample -prefix ${OUTDIR}/tmp.nii -orient ras -input ${OUTDIR}/rev_phase.nii
                  3drefit -deoblique ${OUTDIR}/tmp.nii
                  mv ${OUTDIR}/tmp.nii ${OUTDIR}/rev_phase.nii

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
                  3dTcat -tr ${TR} -prefix ${OUTDIR}/topup_data.nii ${OUTDIR}/blipdown_m.nii[0..4] ${OUTDIR}/blipup_m.nii[0..4]


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
                  #3dAutomask -apply_prefix ${OUTDIR}/tmp.nii ${OUTDIR}/u${PREF}func_data.nii
                  #rm ${OUTDIR}/u${PREF}func_data.nii
                  #mv ${OUTDIR}/tmp.nii ${OUTDIR}/u${PREF}func_data.nii
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

                  gunzip -f ${OUTDIR}/u${PREF}func_data.nii.gz
                  python ${QCDIR}/save_image_diff.py -o ${OUTDIR} -i ${OUTDIR} -a ${PREF}func_data.nii -b u${PREF}func_data.nii -type mean -msg1 'Before FMAP' -msg2 'After FMAP'

                  PREF=u${PREF}
                  log_command_div 'FIELDMAP CORRECTION'
            fi


            # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



            if [ "$DO_ING" -eq "1" ]; then
                  print_debug 'Global intensity normalization' 'fslmaths ${OUTDIR}/${PREF}func_data -ing ${GLOB_VAL} ${OUTDIR}/i${PREF}func_data'
                  fslmaths ${OUTDIR}/${PREF}func_data -thrp 10 -ing ${GLOB_VAL} ${OUTDIR}/i${PREF}func_data
                  gunzip -f ${OUTDIR}/i${PREF}func_data.nii.gz
                  PREF=i${PREF}
            fi


            # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            # Update brain mask
            rm ${OUTDIR}/nat_mask.nii
            #3dAutomask -prefix ${OUTDIR}/nat_mask.nii -apply_prefix ${OUTDIR}/proc_data_native.nii ${OUTDIR}/${PREF}func_data.nii
            3dAutomask -prefix ${OUTDIR}/nat_mask.nii ${OUTDIR}/${PREF}func_data.nii
            cp -f ${OUTDIR}/${PREF}func_data.nii ${OUTDIR}/proc_data_native.nii
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
            fi


      } &> ${OUTDIR}/01_moco.log
      END=$(date -u +%s.%N)
      DIFF=`echo "( $END - $START )" | bc`
      printf "DONE [%.1f s]\n" $DIFF

fi

exit

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
                cp ${OUTDIR}/proc_data_native.nii ${ANATDIR}/func_data.nii

                NVOLS=`fslnvols ${ANATDIR}/func_data.nii`

                #FUNC
                fslmaths ${ANATDIR}/func_data.nii -Tmean -thrp 30 ${ANATDIR}/mean_func_data.nii
                gunzip -f ${ANATDIR}/mean_func_data.nii.gz
                ${ABIN}/DenoiseImage -i ${ANATDIR}/mean_func_data.nii -r 3x3x3 -n Gaussian -o [ ${ANATDIR}/mean_func_data_nds.nii ]

                #ANAT
                3dresample -prefix ${ANATDIR}/anat_orig.nii -orient ras -input ${ANATDIR}/anat.nii
                3drefit -deoblique ${ANATDIR}/anat_orig.nii
                ${ABIN}/N4BiasFieldCorrection -i ${ANATDIR}/anat_orig.nii -s 3  -o [ ${ANATDIR}/anat_wb.nii ]
                python ${WDIR}/extract_brain.py -i "${ANATDIR}/anat_wb.nii" -o "${ANATDIR}/brain_mask.nii"
                ${ABIN}/DenoiseImage -i ${ANATDIR}/anat_wb.nii -p 1x1x1 -r 4x4x4 -n Rician -x ${ANATDIR}/brain_mask.nii -o [ ${ANATDIR}/anat_wbn.nii ]
                3dcalc -a ${ANATDIR}/anat_wbn.nii -b ${ANATDIR}/brain_mask.nii -expr 'a*b' -prefix ${ANATDIR}/anat_proc_brain.nii


                mv ${ANATDIR}/anat_wbn.nii ${ANATDIR}/anat_proc.nii
                copy_header ${ANATDIR}/anat_orig.nii ${ANATDIR}/anat_proc.nii
                rm -f ${ANATDIR}/anat_w*.nii
                rm -f ${ANATDIR}/anat.nii



                copy_header ${ANATDIR}/mean_func_data.nii ${ANATDIR}/mean_func_data_nds.nii

                matlab "-nodesktop -nosplash " <<<"coreg_normalise('${ANATDIR}/mean_func_data_nds.nii', '${ANATDIR}/func_data.nii', ${NVOLS}, '${ANATDIR}/anat_proc_brain.nii', [1 0 0], {'${DARTEL_TEMPLATE_PREF}1.nii', '${DARTEL_TEMPLATE_PREF}2.nii', '${DARTEL_TEMPLATE_PREF}3.nii', '${DARTEL_TEMPLATE_PREF}4.nii', '${DARTEL_TEMPLATE_PREF}5.nii', '${DARTEL_TEMPLATE_PREF}6.nii'}, [5 5 5]); exit;"

                copy_header ${ANATDIR}/anat_proc_brain.nii ${ANATDIR}/anat_proc.nii

                matlab "-nodesktop -nosplash " <<<"coreg_normalise('${ANATDIR}/mean_func_data_nds.nii', '${ANATDIR}/func_data.nii', ${NVOLS}, '${ANATDIR}/anat_proc.nii', [0 1 0], {'${DARTEL_TEMPLATE_PREF}1.nii', '${DARTEL_TEMPLATE_PREF}2.nii', '${DARTEL_TEMPLATE_PREF}3.nii', '${DARTEL_TEMPLATE_PREF}4.nii', '${DARTEL_TEMPLATE_PREF}5.nii', '${DARTEL_TEMPLATE_PREF}6.nii'}, [5 5 5]); exit;"

                #matlab "-nodesktop -nosplash " <<<"affine_spm2fsl('${ANATDIR}', '${ANATDIR}/anat_proc.nii','${ANATDIR}/mean_func_data_nds.nii'); exit;"

          fi
          rm -f ${ANATDIR}/anat_native.nii

          mv ${ANATDIR}/rc1anat_proc.nii ${ANATDIR}/rc1anat_proc_dartel.nii
          mv ${ANATDIR}/rc2anat_proc.nii ${ANATDIR}/rc2anat_proc_dartel.nii
          mv ${ANATDIR}/rc3anat_proc.nii ${ANATDIR}/rc3anat_proc_dartel.nii

          matlab "-nodesktop -nosplash " <<<"reslice_images('${ANATDIR}/mean_func_data_nds.nii', {'${ANATDIR}/anat_proc.nii', '${ANATDIR}/anat_proc_brain.nii', '${ANATDIR}/c1anat_proc.nii', '${ANATDIR}/c2anat_proc.nii', '${ANATDIR}/c3anat_proc.nii'}); exit;"


          mv ${ANATDIR}/ranat_proc.nii ${ANATDIR}/anat_native.nii
          mv ${ANATDIR}/ranat_proc_brain.nii ${ANATDIR}/anat_native_brain.nii

          mv ${ANATDIR}/rc1anat_proc.nii ${ANATDIR}/gm_tpm_native.nii
          mv ${ANATDIR}/rc2anat_proc.nii ${ANATDIR}/wm_tpm_native.nii
          mv ${ANATDIR}/rc3anat_proc.nii ${ANATDIR}/csf_tpm_native.nii


          copy_header ${ANATDIR}/mean_func_data_nds.nii ${ANATDIR}/anat_native.nii
          copy_header ${ANATDIR}/mean_func_data_nds.nii ${ANATDIR}/anat_native_brain.nii
          copy_header ${ANATDIR}/mean_func_data_nds.nii ${ANATDIR}/gm_tpm_native.nii
          copy_header ${ANATDIR}/mean_func_data_nds.nii ${ANATDIR}/wm_tpm_native.nii
          copy_header ${ANATDIR}/mean_func_data_nds.nii ${ANATDIR}/csf_tpm_native.nii



          flirt -in ${ANATDIR}/anat_native.nii -ref ${ANATDIR}/anat_proc.nii -omat ${ANATDIR}/func2anat_fsl.mat -out ${ANATDIR}/anat_native_fsl \
               -usesqform -dof 6 -coarsesearch 18 -finesearch 9 -searchrx -20 20 -searchry -20 20 -searchrz -20 20


          convert_xfm -omat ${ANATDIR}/anat2func_fsl.mat -inverse ${ANATDIR}/func2anat_fsl.mat
          flirt -applyxfm -init ${ANATDIR}/anat2func_fsl.mat -in ${ANATDIR}/anat_proc.nii -out ${ANATDIR}/anat_native_fsl_func.nii.gz -ref ${ANATDIR}/mean_func_data.nii


      if [ "$DO_CLEAN" -eq "1" ]; then

            rm -f ${ANATDIR}/c1anat_proc.nii
            rm -f ${ANATDIR}/c2anat_proc.nii
            rm -f ${ANATDIR}/c3anat_proc.nii
            rm -f ${ANATDIR}/c4anat_proc.nii
            rm -f ${ANATDIR}/c5anat_proc.nii
            rm -f ${ANATDIR}/c6anat_proc.nii

            rm -f ${ANATDIR}/anat.nii
            rm -f ${ANATDIR}/func_data.nii
            rm -f ${ANATDIR}/func_data.mat
            rm -f ${ANATDIR}/anat_wb.nii
            rm -f ${ANATDIR}/brain_mask.nii
            gzip -f ${ANATDIR}/iy_anat_proc.nii
            gzip -f ${ANATDIR}/y_anat_proc.nii

      fi
      } &> ${OUTDIR}/02_ant_norm.log

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
            3dAutomask -prefix ${OUTDIR}/nat_mask.nii -peel 2 ${OUTDIR}/proc_data_native.nii
            3dTproject -prefix ${OUTDIR}/tmp_melodic.nii -polort 1 -stopband 0 0.01 -mask ${OUTDIR}/nat_mask.nii  -TR ${TR} -input ${OUTDIR}/proc_data_native.nii -blur 6

            rm -f -R ${OUTDIR}/melodic.ic
            rm -f -R ${OUTDIR}/FIX

            copy_header ${OUTDIR}/proc_data_native.nii ${OUTDIR}/tmp_melodic.nii

            3dpc -eigonly -vmean -vnorm -mask ${OUTDIR}/nat_mask.nii -prefix ${OUTDIR}/pc_var ${OUTDIR}/tmp_melodic.nii
            IDX=0
            NIC=1
            PREV=-1
            STABILITY=0
            while IFS=" " read -r value1 value2 value3 value4
            do
                  if (( IDX > 0 )); then
                        if (( $(echo "${value4} > 0.9" | bc -l) )); then
                          PCT=$(echo "${value4} * 100" | bc -l)
                          NIC=${IDX}
                          break
                        fi

                        if (( $(echo "${PREV} < 0" | bc -l) )); then
                              PREV=${value4}
                              PREV_CHANGE=0
                        else

                              CHG=$(echo "${value4} - ${PREV}" | bc -l)
                              UPD=$(echo "${CHG} - ${PREV_CHANGE}" | bc -l)

                              if (( $(echo "${UPD} < 0" | bc -l) )); then
                                    UPD=$(echo "${UPD}*-1" | bc -l )
                              fi

                              if (( $(echo "${UPD} < 0.00001" | bc -l) )); then
                                echo "[$IDX] Change = ${UPD}*"
                              else
                                echo "[$IDX] Change = ${UPD}"
                              fi
                              if (( $(echo "${UPD} < 0.00001" | bc -l) )); then
                                  if (( $STABILITY > 5 )); then
                                        echo "[$IDX] Change = ${UPD}"
                                        PCT=$(echo "${value4} * 100" | bc -l)
                                        NIC=$IDX
                                        break
                                  else
                                        STABILITY=$((STABILITY+1))
                                  fi
                              fi

                              PREV=${value4}
                              PREV_CHANGE=${CHG}
                        fi

                  fi

                  IDX=$((IDX+1))
            done < "${OUTDIR}/pc_var_eig.1D"

            echo "EXTRACTING ${NIC} components (${PCT}% of variance explained) using ${ICA_TYPE} estimation"

            # --dimest=lap is the default, but it seems empirically to overestimate components
            # see also Varoquaux et al. 2010 NeuroImage
            if [ "$ICA_TYPE" == "melodic" ]; then
                  melodic -i ${OUTDIR}/tmp_melodic.nii -o ${OUTDIR}/melodic.ic --tr=${TR} --mmthresh=0.5 --nobet --dim=${NIC} --nl=pow4 \
                         --mask=${OUTDIR}/nat_mask.nii --Ostats --report --eps=0.0001 --bgimage=${ANATDIR}/anat_native_fsl_func
            fi

            if [ "$ICA_TYPE" == "canica" ]; then
                  mkdir ${OUTDIR}/melodic.ic
                  #mkdir ${OUTDIR}/melodic.ic/report

                  #if test -f "${OUTDIR}/Group_Comps_Native.nii"; then
                  #      python ${WDIR}/canica_estimation.py -o ${OUTDIR}/melodic.ic -in ${OUTDIR}/tmp_melodic.nii \
                  #              -mask ${OUTDIR}/nat_mask.nii -nIC ${NIC} -decomp dictlearning -dict_init ${OUTDIR}/Group_Comps_Native.nii
                  #else
                        python ${WDIR}/canica_estimation.py -o ${OUTDIR}/melodic.ic -in ${OUTDIR}/tmp_melodic.nii \
                                -mask ${OUTDIR}/nat_mask.nii -nIC ${NIC} #-decomp dictlearning
                  #fi

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

            fi


            rm ${OUTDIR}/tmp_melodic.nii

            # Creates a smoothed copy of the melodic components. For visualization purposes only
            3dmerge -doall -1blur_fwhm 5 -prefix ${OUTDIR}/melodic.ic/melodic_IC_smooth.nii ${OUTDIR}/melodic.ic/melodic_IC.nii.gz


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
            cp ${OUTDIR}/nat_mask.nii ${FIXDIR}/mask.nii
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
      } &> ${OUTDIR}/04_ICA.log

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

      SUFFIX=''
      if test -f "${OUTDIR}/proc_data_native_nlm.nii"; then
            SUFFIX='_nlm'
      fi
      copy_header ${ANATDIR}/mean_func_norm.nii  ${ANATDIR}/csf_tpm_native.nii
      copy_header ${ANATDIR}/mean_func_norm.nii  ${ANATDIR}/gm_tpm_native.nii
      copy_header ${ANATDIR}/mean_func_norm.nii  ${ANATDIR}/wm_tpm_native.nii

      # ERODE tissue masks
      fslmaths ${ANATDIR}/csf_tpm_native.nii -kernel 2D -fmedian -thr 0.8 -bin ${ANATDIR}/csf_mask_native.nii
      gunzip -f ${ANATDIR}/csf_mask_native.nii.gz

      fslmaths ${ANATDIR}/wm_tpm_native.nii -kernel 2D -fmedian -thr 0.8 -bin ${ANATDIR}/wm_mask_native.nii
      gunzip -f ${ANATDIR}/wm_mask_native.nii.gz

      fslmaths ${ANATDIR}/gm_tpm_native.nii -kernel 2D -fmedian -thr 0.5 -bin ${ANATDIR}/gm_mask_native.nii
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

      } &> ${OUTDIR}/03_SignalExtraction.log

      END=$(date -u +%s.%N)
      DIFF=`echo "( $END - $START )" | bc`
      printf "DONE [%.1f s]\n" $DIFF

fi

# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------


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
        source ${FIXBIN}/fix -a ${FIXDIR}/fix4melview_${FIX_CL_LABEL}_thr${FIX_THR}.txt
        cd ${CURDIR}

        mv ${FIXDIR}/filtered_func_data_clean.nii.gz ${OUTDIR}/proc_data_native_fix.nii.gz
        gunzip -f ${OUTDIR}/proc_data_native_fix.nii.gz
        python ${QCDIR}/save_image_diff.py -o ${OUTDIR} -i ${OUTDIR} -a proc_data_native${SUFFIX}.nii -b  proc_data_native_fix.nii -type std -msg1 'Before FIX' -msg2 'After FIX'

        mv ${OUTDIR}/*.png ${IMGDIR}
        mv ${OUTDIR}/*.gif ${IMGDIR}
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

        rm -f -r ${OUTDIR}/dicer
        rm -f -r ${OUTDIR}/dicer_fix
        CURDIR=`pwd`

        # TODO Set up conf for DiCER PATH
        cd ${DICERPATH}


        mkdir ${OUTDIR}/dicer_fix
        #mkdir ${OUTDIR}/dicer

        cp ${OUTDIR}/motion_regressors_12_24.txt ${OUTDIR}/dicer_fix/regressors.txt
        cp ${OUTDIR}/motion_regressors_12_24.txt ${OUTDIR}/dicerregressors.txt

        # If background is not 0, fast has problem getting the correct tissue distrub
        #3dcalc -a ${OUTDIR}/proc_data_native_fix.nii -b ${OUTDIR}/gm_mask_nat.nii -expr "a*b" -prefix ${OUTDIR}/tmp.nii
        cp -f ${OUTDIR}/proc_data_native_fix.nii ${OUTDIR}/tmp.nii

        fslmaths ${OUTDIR}/tmp.nii -Tmean ${OUTDIR}/tmp_mean.nii
        gunzip -f ${OUTDIR}/tmp_mean.nii.gz


        flirt -in ${ANATDIR}/anat_native_brain.nii -ref ${ANATDIR}/anat_proc.nii  -out ${ANATDIR}/anat_native_brain_fsl \
             -usesqform -dof 6 -coarsesearch 18 -finesearch 9 -searchrx -20 20 -searchry -20 20 -searchrz -20 20
        flirt -applyxfm -init ${ANATDIR}/anat2func_fsl.mat -in ${ANATDIR}/anat_proc_brain.nii -out ${ANATDIR}/anat_native_brain_fsl_func.nii.gz -ref ${ANATDIR}/mean_func_data.nii
        gunzip -f ${ANATDIR}/anat_native_brain_fsl_func.nii.gz
        cp -f ${ANATDIR}/anat_native_brain_fsl_func.nii ${ANATDIR}/tmp_anat_native.nii
        copy_header ${OUTDIR}/tmp_mean.nii ${ANATDIR}/tmp_anat_native.nii
        rm -f ${OUTDIR}/tmp_mean.nii

        gzip -f ${OUTDIR}/tmp.nii
        rm -f ${OUTDIR}/tmp.nii

        source DiCER_lightweight.sh -i ${OUTDIR}/tmp.nii.gz -a ${ANATDIR}/tmp_anat_native.nii -w ${OUTDIR}/dicer_fix -s SUBJECT_1_FIX  -p 1 -d

        cd ${CURDIR}
        rm ${ANATDIR}/tmp_anat_native.nii
        #mv ${OUTDIR}/dicer_fix/tmp_detrended_hpf_dbscan.nii.gz  ${OUTDIR}/proc_data_native_fix_d.nii.gz
        #gunzip -f ${OUTDIR}/proc_data_native_fix_d.nii.gz

        cp -f ${OUTDIR}/dicer_fix/SUBJECT_1_FIX_dbscan_liberal_regressors.tsv ${OUTDIR}/DiCER_regressors.txt

        rm -f -r ${OUTDIR}/dicer_fix
        rm -f ${OUTDIR}/tmp.nii.gz
        rm -f ${OUTDIR}/tmp.nii
        rm -f ${OUTDIR}/proc_data_native_fix_d.nii
        #python ${QCDIR}/save_image_diff.py -o ${OUTDIR} -i ${OUTDIR} -a proc_data_native_fix.nii -b  proc_data_native_fix_d.nii -type std -msg1 'Before DiCER' -msg2 'After DiCER'
        #mv ${OUTDIR}/*.png ${IMGDIR}
        #mv ${OUTDIR}/*.gif ${IMGDIR}

      }  &> ${OUTDIR}/08_DiCER.log
      END=$(date -u +%s.%N)
      DIFF=`echo "( $END - $START )" | bc`
      printf "DONE [%.1f s]\n" $DIFF
fi

# Apply normalisation WARP
if [ "$DO_NORM" -eq "1" ]; then

      printf "[$SUB] Warping functional images to MNI space ... "
      START=$(date -u +%s.%N)

      {

            rm -f ${ANATDIR}/c6anat_proc.nii
            # Run selected nuisance regression
            SUF=(mni_fix mni)
            IDX=0

            rm -f ${OUTDIR}/proc_data_mni*.nii


            for INPREF in proc_data_native_fix proc_data_native
            do
                CURR_SUFF=${SUF[$IDX]}
                IDX=$((IDX+1))
                if ! test -f "${OUTDIR}/${INPREF}.nii"; then
                      continue
                fi

                rm -f ${OUTDIR}/proc_data_native_fix_tmp.nii
                rm -f ${OUTDIR}/rproc_data_native_fix_tmp.nii
                rm -f ${OUTDIR}/rproc_data_native_fix_tmp.mat
                rm -f ${OUTDIR}/proc_data_native_fix_tmp.mat
                rm -f ${OUTDIR}/mean_func_data_nds.nii

                rm -f ${OUTDIR}/fix_mean.nii

                cp -f ${OUTDIR}/${INPREF}.nii ${OUTDIR}/proc_data_native_fix_tmp.nii
                NVOLS=`fslnvols ${OUTDIR}/proc_data_native_fix_tmp.nii`
                3dTcat -prefix ${OUTDIR}/fix_mean.nii ${OUTDIR}/proc_data_native_fix_tmp.nii[0]
                cp -f ${ANATDIR}/mean_func_data_nds.nii ${OUTDIR}/mean_func_data_nds.nii

                matlab "-nodesktop -nosplash " <<<"coreg_same_image('${OUTDIR}/mean_func_data_nds.nii', '${OUTDIR}/fix_mean.nii', '${OUTDIR}/proc_data_native_fix_tmp.nii', ${NVOLS}); exit;"
                matlab "-nodesktop -nosplash " <<<"coreg_normalise('${OUTDIR}/mean_func_data_nds.nii', '${OUTDIR}/proc_data_native_fix_tmp.nii', ${NVOLS}, '${ANATDIR}/anat_proc.nii', [0 0 1], {'${DARTEL_TEMPLATE_PREF}1.nii', '${DARTEL_TEMPLATE_PREF}2.nii', '${DARTEL_TEMPLATE_PREF}3.nii', '${DARTEL_TEMPLATE_PREF}4.nii', '${DARTEL_TEMPLATE_PREF}5.nii', '${DARTEL_TEMPLATE_PREF}6.nii'}, [0 0 0]); exit;"

                cp ${ANATDIR}/rc3anat_proc.nii ${ANATDIR}/csf_tpm_native_.nii
                cp ${ANATDIR}/rc2anat_proc.nii ${ANATDIR}/wm_tpm_native_.nii
                cp ${ANATDIR}/rc1anat_proc.nii ${ANATDIR}/gm_tpm_native_.nii

                matlab "-nodesktop -nosplash " <<<"coreg_normalise('${OUTDIR}/mean_func_data_nds.nii', '${ANATDIR}/csf_tpm_native.nii', ${NVOLS}, '${ANATDIR}/anat_proc.nii', [0 0 1], {'${DARTEL_TEMPLATE_PREF}1.nii', '${DARTEL_TEMPLATE_PREF}2.nii', '${DARTEL_TEMPLATE_PREF}3.nii', '${DARTEL_TEMPLATE_PREF}4.nii', '${DARTEL_TEMPLATE_PREF}5.nii', '${DARTEL_TEMPLATE_PREF}6.nii'}, [0 0 0]); exit;"
                matlab "-nodesktop -nosplash " <<<"coreg_normalise('${OUTDIR}/mean_func_data_nds.nii', '${ANATDIR}/wm_tpm_native.nii', ${NVOLS}, '${ANATDIR}/anat_proc.nii', [0 0 1], {'${DARTEL_TEMPLATE_PREF}1.nii', '${DARTEL_TEMPLATE_PREF}2.nii', '${DARTEL_TEMPLATE_PREF}3.nii', '${DARTEL_TEMPLATE_PREF}4.nii', '${DARTEL_TEMPLATE_PREF}5.nii', '${DARTEL_TEMPLATE_PREF}6.nii'}, [0 0 0]); exit;"
                matlab "-nodesktop -nosplash " <<<"coreg_normalise('${OUTDIR}/mean_func_data_nds.nii', '${ANATDIR}/gm_tpm_native.nii', ${NVOLS}, '${ANATDIR}/anat_proc.nii', [0 0 1], {'${DARTEL_TEMPLATE_PREF}1.nii', '${DARTEL_TEMPLATE_PREF}2.nii', '${DARTEL_TEMPLATE_PREF}3.nii', '${DARTEL_TEMPLATE_PREF}4.nii', '${DARTEL_TEMPLATE_PREF}5.nii', '${DARTEL_TEMPLATE_PREF}6.nii'}, [0 0 0]); exit;"

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
                mv ${OUTDIR}/wproc_data_native_fix_tmp.nii ${OUTDIR}/proc_data_mni.nii
                mv ${ANATDIR}/wcsf_tpm_native.nii ${ANATDIR}/csf_mask_group.nii
                mv ${ANATDIR}/wwm_tpm_native.nii ${ANATDIR}/wm_mask_group.nii
                mv ${ANATDIR}/wgm_tpm_native.nii ${ANATDIR}/gm_mask_group.nii


                fslmaths ${ANATDIR}/csf_mask_group.nii -thr 0.3 -bin -ero  ${ANATDIR}/csf_mask_group_.nii
                gunzip -f ${ANATDIR}/csf_mask_group_.nii.gz
                mv ${ANATDIR}/csf_mask_group_.nii  ${ANATDIR}/csf_mask_group.nii


                fslmaths ${ANATDIR}/wm_mask_group.nii -thr 0.4 -bin  -ero  ${ANATDIR}/wm_mask_group_.nii
                gunzip -f ${ANATDIR}/wm_mask_group_.nii.gz
                mv ${ANATDIR}/wm_mask_group_.nii  ${ANATDIR}/wm_mask_group.nii

                fslmaths ${ANATDIR}/gm_mask_group.nii -thr 0.5 -bin   ${ANATDIR}/gm_mask_group_.nii
                gunzip -f ${ANATDIR}/gm_mask_group_.nii.gz
                mv ${ANATDIR}/gm_mask_group_.nii  ${ANATDIR}/gm_mask_group.nii

                3dAutomask -prefix ${OUTDIR}/mask_mni.nii ${OUTDIR}/proc_data_mni.nii
                3dBlurToFWHM -quiet -prefix ${OUTDIR}/proc_data_mni_s.nii -input ${OUTDIR}/proc_data_mni.nii -FWHM 5 -rate 2 -mask ${OUTDIR}/mask_mni.nii

                mv ${OUTDIR}/proc_data_mni_s.nii ${OUTDIR}/proc_data_mni.nii

                fslmaths  ${OUTDIR}/proc_data_mni.nii -ing 1000 ${OUTDIR}/proc_data_mni.nii -odt short
                gunzip -f ${OUTDIR}/proc_data_mni.nii.gz
                mv ${OUTDIR}/proc_data_mni.nii ${OUTDIR}/proc_data_${CURR_SUFF}.nii
            done


} &> ${OUTDIR}/05_norm2mni.log

      END=$(date -u +%s.%N)
      DIFF=`echo "( $END - $START )" | bc`
      printf "DONE [%.1f s]\n" $DIFF
fi

if [ "$DO_ATLAS_NAT" -eq "1" ]; then
    printf "[$SUB] MNI Atlas -> Subject space ... "
    START=$(date -u +%s.%N)
    {
        #
        WDIR='/home/luna.kuleuven.be/u0101486/workspace/fmri_proc/'
        for atlas in aal cobra hammers local_global_400 neuromorphometrics lg400_cobra; do
            echo $atlas
            cp -f ${WDIR}/atlas/CAT12/${atlas}.txt ${ANATDIR}/${atlas}.txt
            ${ABIN}/antsApplyTransforms \
                  -i ${WDIR}/atlas/CAT12/${atlas}.nii \
                  --float \
                  -r ${DARTEL_TEMPLATE_DIR}/Group_T1_Avg_Brain.nii \
                  -o ${ANATDIR}/${atlas}_group.nii \
                  -t [${DARTEL_TEMPLATE_DIR}/mni2group0GenericAffine.mat,0] \
                  -t [${DARTEL_TEMPLATE_DIR}/mni2group1Warp.nii.gz,0] \
                  -n NearestNeighbor

            matlab "-nodesktop -nosplash " <<<"group2sub( '${ANATDIR}/u_rc1anat_proc.nii', '${ANATDIR}/mean_func_data_nds.nii',{'${ANATDIR}/${atlas}_group.nii'}, '${ANATDIR}/w${atlas}_group_u_rc1anat_proc.nii'); exit;"
            rm -f ${ANATDIR}/w${atlas}_group_u_rc1anat_proc.nii
            mv  ${ANATDIR}/rw${atlas}_group_u_rc1anat_proc.nii ${ANATDIR}/${atlas}_native.nii

            cp ${ANATDIR}/mean_func_norm.nii ${ANATDIR}/mean_func_norm_tmp.nii
            cp ${ANATDIR}/mean_func_data_nds.nii ${ANATDIR}/mean_func_data_tmp.nii
            matlab "-nodesktop -nosplash " <<<"coreg_same_image('${ANATDIR}/mean_func_norm_tmp.nii', '${ANATDIR}/mean_func_data_tmp.nii', {'${ANATDIR}/${atlas}_native.nii'}, 0); exit;"
            rm -f ${ANATDIR}/mean_func_norm_tmp.nii
            rm -f ${ANATDIR}/mean_func_data_tmp.nii
            rm -f ${ANATDIR}/rmean_func_data_tmp.nii
            mv ${ANATDIR}/r${atlas}_native.nii ${ANATDIR}/${atlas}_native.nii
            rm -f ${ANATDIR}/${atlas}_group.nii
        done
    }  &> ${OUTDIR}/09_Atlas_norm.log
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

if [ "$DO_QA" -eq "1" ]; then
      # Performs QS in native space
      printf "[$SUB] Generating QC plots ... "
      START=$(date -u +%s.%N)

      {
          #fslmaths ${OUTDIR}/proc_data_native_nlm.nii -Tmean ${OUTDIR}/func_mean.nii
          #gunzip -f ${OUTDIR}/func_mean.nii.gz

          cp -f ${DARTEL_TEMPLATE_DIR}/Group_Mask.nii ${ANATDIR}/GroupMask.nii


          model_qa()
          {
              OUTDIR=${1}
              ANATDIR=${2}
              MNI_REF=${3}
              QAOUT=${4}
              TR=${5}
              QCDIR=${6}
              INPREF=${7}
              REG_MODEL=${8}

              mkdir -p ${OUTDIR}/${QAOUT}/QA_${REG_MODEL}
              QA_DIR=${OUTDIR}/${QAOUT}/QA_${REG_MODEL}

              rm -f ${OUTDIR}/${INPREF}_${REG_MODEL}.nii
              rm -f ${OUTDIR}/${INPREF}_${REG_MODEL}.nii.gz
              rm -f ${QA_DIR}/*


              fslmeants -i ${OUTDIR}/${INPREF}.nii -m ${ANATDIR}/csf_mask_group.nii -o ${OUTDIR}/csf_sig_${REG_MODEL}.txt
              fslmeants -i ${OUTDIR}/${INPREF}.nii -m ${ANATDIR}/wm_mask_group.nii -o ${OUTDIR}/wm_sig_${REG_MODEL}.txt
              3dcalc -a ${ANATDIR}/csf_mask_group.nii -b ${ANATDIR}/wm_mask_group.nii -expr 'a+b' -prefix ${ANATDIR}/nongm_mask_group_${REG_MODEL}.nii
              CMD="3dTproject -input ${OUTDIR}/${INPREF}.nii -TR ${TR} -norm  "
              CMD="${CMD}  "
              CMD="${CMD} -prefix ${OUTDIR}/${INPREF}_${REG_MODEL}.nii "
              CMD="${CMD} -cenmode NTRP -censor ${OUTDIR}/temporal_mask_fd.txt "


              if [ "${REG_MODEL}" == "SRP24WM1CSF1" ] || [ "${REG_MODEL}" == "SRP24WM1CSF1_D" ] || \
                 [ "${REG_MODEL}" == "SRP24CC" ]; then
                  CMD="${CMD} -ort ${OUTDIR}/motion_regressors_12.txt "
              fi

              if [ "${REG_MODEL}" == "SFIX_CC" ]; then
                  rm -f ${OUTDIR}/dvars.txt

                  python run_acompcor.py -d ${OUTDIR} -i ${OUTDIR}/${INPREF}.nii -n ${ANATDIR}/nongm_mask_group_${REG_MODEL}.nii -b ${OUTDIR}/mask_mni.nii -t ${TR} \
                                        -aout acompcor_fix -tout tcompcor_fix -var 0.6

                  # Remove headers from compcor results
                  sed '1d' ${OUTDIR}/acompcor_fix.txt > ${OUTDIR}/tmp_acompcor_fix.txt
                  mv ${OUTDIR}/tmp_acompcor_fix.txt ${OUTDIR}/acompcor_fix.txt

                  sed '1d' ${OUTDIR}/tcompcor_fix.txt > ${OUTDIR}/tmp_tcompcor_fix.txt
                  mv ${OUTDIR}/tmp_tcompcor_fix.txt ${OUTDIR}/tcompcor_fix.txt

                  python fix_acompcor_pc.py -compcor ${OUTDIR}/acompcor_fix.txt \
                      -fix_ics ${OUTDIR}/tcompcor_fix.txt -expVar 0.8 -svar 0.01 \
                      -mot12 ${OUTDIR}/motion_regressors_12.txt \
                      -mot24 ${OUTDIR}/motion_regressors_24.txt \
                      -out ${OUTDIR}/fix_cc_nuis.txt

                  #1dBandpass -norm -dt ${TR} 0.009 100 ${OUTDIR}/fix_cc_nuis.txt > ${OUTDIR}/fix_cc_nuis_bp.txt
                  CMD="${CMD} -ort ${OUTDIR}/ic_tcs.txt "
                  CMD="${CMD} -ort ${OUTDIR}/acompcor_fix.txt "
                  #CMD="${CMD} -ort ${OUTDIR}/fix_cc_nuis.txt "
                  #CMD="${CMD} -ort ${OUTDIR}/dvars.txt "
              fi

              if [ "${REG_MODEL}" == "SFIX" ] ; then
                  CMD="${CMD} -ort ${OUTDIR}/csf_sig_${REG_MODEL}.txt "
                  CMD="${CMD} -ort ${OUTDIR}/wm_sig_${REG_MODEL}.txt "
              fi

              if [ "${REG_MODEL}" == "SFIX_D" ] || [ "${REG_MODEL}" == "SRP24WM1CSF1_D" ]; then
                  CMD="${CMD} -ort ${OUTDIR}/DiCER_regressors.txt "
              fi

              if [ "${REG_MODEL}" == "SRP9" ] || [ "${REG_MODEL}" == "SRP24WM1CSF1" ]; then
                  CMD="${CMD} -ort ${OUTDIR}/csf_sig_${REG_MODEL}.txt "
                  CMD="${CMD} -ort ${OUTDIR}/wm_sig_${REG_MODEL}.txt "
              fi



              if [ "${REG_MODEL}" == "SRP24WM1CSF1" ] || [ "${REG_MODEL}" == "SRP24CC" ] || \
                 [ "${REG_MODEL}" == "SRP24WM1CSF1_D" ]; then
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
                  CMD="${CMD} -passband 0.01 0.15 -polort 1 "
              else
                  CMD="${CMD} -passband 0.01 0.15 -polort 0 "
              fi

              echo "RUNNING: ${CMD}"
              eval "${CMD}"
              rm -f ${OUTDIR}/csf_${REG_MODEL}.txt
              rm -f ${OUTDIR}/wm_${REG_MODEL}.txt


              #NVOLS=`fslnvols ${OUTDIR}/${INPREF}_${REG_MODEL}.nii`
              #matlab "-nodesktop -nosplash " <<<"coreg_normalise('${ANATDIR}/mean_func_norm.nii', '${OUTDIR}/${INPREF}_${REG_MODEL}.nii', ${NVOLS}, '${ANATDIR}/anat_proc.nii', [0 0 1], {'${DARTEL_TEMPLATE_PREF}1.nii', '${DARTEL_TEMPLATE_PREF}2.nii', '${DARTEL_TEMPLATE_PREF}3.nii', '${DARTEL_TEMPLATE_PREF}4.nii', '${DARTEL_TEMPLATE_PREF}5.nii', '${DARTEL_TEMPLATE_PREF}6.nii'}, [5 5 5]); exit;"
              3dReHo  -inset ${OUTDIR}/${INPREF}_${REG_MODEL}.nii  -prefix ${QA_DIR}/REHO.nii
              #rm -f ${ANATDIR}/w${INPREF}_${REG_MODEL}.nii

              python ${QCDIR}/QC_grey_plot.py -out ${QA_DIR} -in ${OUTDIR}  -mpe motion_estimate.par \
                                -fname ${INPREF}_${REG_MODEL}.nii -csf_name ${ANATDIR}/csf_mask_group.nii \
                                -wm_name ${ANATDIR}/wm_mask_group.nii \
                                -gm_name ${ANATDIR}/gm_mask_group.nii \
                                -norm zscore -range 1.0 \
                                -prog AFNI -outf 02_Greyplot_${REG_MODEL}  -dpi 100 -tr ${TR}

          #    python ${QCDIR}/QC_FC.py -out ${QA_DIR} -in ${OUTDIR}  \
            #              -fname ${INPREF}_${REG_MODEL}.nii -outf FC_${REG_MODEL}  \
            #              -atlas ${ANATDIR}/aal_native.nii \
            #              -atlas_name aal \
            #              -vmin -0.6 -vmax 0.6 \
            #               -dpi 72

             #python ${QCDIR}/QC_FC.py -out ${QA_DIR} -in ${OUTDIR}  \
            #             -fname ${INPREF}_${REG_MODEL}.nii -outf FC_${REG_MODEL}  \
            #             -atlas ${ANATDIR}/neuromorphometrics_native.nii \
            #             -atlas_name neuromorphometrics \
            #             -vmin -0.6 -vmax 0.6 \
            #              -dpi 72


            #  python ${QCDIR}/QC_FC.py -out ${QA_DIR} -in ${OUTDIR}  \
            #              -fname ${INPREF}_${REG_MODEL}.nii -outf FC_${REG_MODEL}  \
            #              -atlas ${ANATDIR}/cobra_native.nii \
            #              -atlas_name cobra \
            #              -vmin -0.6 -vmax 0.6 \
            #               -dpi 72

            # python ${QCDIR}/QC_FC.py -out ${QA_DIR} -in ${OUTDIR}  \
            #             -fname ${INPREF}_${REG_MODEL}.nii -outf FC_${REG_MODEL}  \
            #             -atlas ${ANATDIR}/hammers_native.nii \
            #             -atlas_name hammers \
            #             -vmin -0.6 -vmax 0.6 \
            #              -dpi 72

              python ${QCDIR}/QC_FC.py -out ${QA_DIR} -in ${OUTDIR}  \
                          -fname ${INPREF}_${REG_MODEL}.nii -outf FC_${REG_MODEL}  \
                          -atlas /home/luna.kuleuven.be/u0101486/workspace/fmri_proc/atlas/CAT12/lg400_cobra.nii \
                          -atlas_name local_global_cobra \
                          -vmin -0.6 -vmax 0.6 \
                           -dpi 72

              rm -f ${OUTDIR}/${INPREF}_${REG_MODEL}.nii
              rm -f ${OUTDIR}/${INPREF}_${REG_MODEL}.nii.gz

              rm -f ${OUTDIR}/sw${INPREF}_${REG_MODEL}.nii
              rm -f ${OUTDIR}/sw${INPREF}_${REG_MODEL}.nii.gz
      }
      export -f model_qa

      parallel -j4 --line-buffer model_qa ::: ${OUTDIR} ::: ${ANATDIR} ::: ${DARTEL_TEMPLATE_PREF} ::: QA_NOFIX ::: ${TR} ::: ${QCDIR}  :::  proc_data_mni ::: SFIX_CC SFIX_D SRP24WM1CSF1 SRP24WM1CSF1_D SRP24CC SRP9
      parallel -j4 --line-buffer model_qa ::: ${OUTDIR} ::: ${ANATDIR} ::: ${DARTEL_TEMPLATE_PREF} ::: QA_FIX ::: ${TR} ::: ${QCDIR}  :::  proc_data_mni_fix ::: SFIX SFIX_D SRP24WM1CSF1 SRP24WM1CSF1_D SRP24CC SRP9
} &> ${OUTDIR}/06_qa.log

      END=$(date -u +%s.%N)
      DIFF=`echo "( $END - $START )" | bc`
      printf "DONE [%.1f s]\n" $DIFF
fi



rm -f ${OUTDIR}/proc_data_native*.nii
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
