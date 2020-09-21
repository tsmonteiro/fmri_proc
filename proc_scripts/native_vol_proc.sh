#!/bin/bash
# Clear terminal screen
#printf "\033c"

OUTDIR=$1
CONFIG=$2

WDIR='/home/luna.kuleuven.be/u0101486/workspace/fmri_proc/'
QCDIR="${WDIR}/QualityControl/"

# Read relevant parameters from CONFIG file
ABIN=$(awk -F\=  '/^ABIN/{print $2}' "${CONFIG}" | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
MOCO=$(awk -F\=  '/^\<MOCO\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')

DO_CLEAN=$(awk -F\=  '/^\<DO_CLEAN\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')

DO_FUNC_BIAS=$(awk -F\=  '/^\<DO_FUNC_BIAS\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')

TR=$(awk -F\=  '/^\<TR\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
ACQ_TYPE=$(awk -F\=  '/^\<ACQ_TYPE\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')

DO_SLC=$(awk -F\=  '/^\<DO_SLC\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
DO_FMAP=$(awk -F\=  '/^\<DO_FMAP\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
EES=$(awk -F\=  '/^\<EES\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
DO_PEST=$(awk -F\=  '/^\<DO_PEST\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')

DO_DEOBL=$(awk -F\=  '/^\<DO_DEOBL\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
DO_DPK=$(awk -F\=  '/^\<DO_DPK\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
DO_ING=$(awk -F\=  '/^\<DO_ING\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
GLOB_VAL=$(awk -F\=  '/^\<GLOB_VAL\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')


DO_FWHM=$(awk -F\=  '/^\<DO_FWHM\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
DO_NLM=$(awk -F\=  '/^\<DO_NLM\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')




print_time()
{
      HOUR=`date +%T`
      DAY=`date +%F`
      TYPE=$1


      echo "[$TYPE] (${DAY} - ${HOUR} )"

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

if [ ! -d "${OUTDIR}/logs" ]; then
      mkdir ${OUTDIR}/logs
fi


NVOLS=`fslnvols ${OUTDIR}/func_data.nii`
NZ=`3dinfo -nk ${OUTDIR}/func_data.nii`

#print_debug "NVOLS = ${NVOLS}"


PREF=''

# Slomoco is better performed, it seems, prior to anything else
if [ "$MOCO" == "slomoco" ]; then
      print_time 'INFO'
      echo "Performing motion correction using SLOMOCO method [${OUTDIR}/logs/moco.log]"
      {
          source ./align_pestica.sh ${OUTDIR} ${PREF}func_data ${TR} 1

          # This is here only for comparison purposes
          3dvolreg -heptic -prefix ${OUTDIR}/r2${PREF}func_data.nii -base 0 -rot_thresh 0.01 -delta 2 \
                -x_thresh 0.01 -zpad 10 -maxite 60 -1Dfile ${OUTDIR}/motion_estimate.par -maxdisp1D ${OUTDIR}/maximum_disp.1d \
                 ${OUTDIR}/${PREF}func_data.nii


          3drefit -deoblique ${OUTDIR}/r2${PREF}func_data.nii
          1dplot -thick -volreg -png ${OUTDIR}/motion_estimate.png -one ${OUTDIR}/motion_estimate.par
          1dplot -thick -png ${OUTDIR}/maximum_disp.png -one ${OUTDIR}/maximum_disp.1d
          1dplot -thick -png ${OUTDIR}/maximum_disp_delt.png -one ${OUTDIR}/maximum_disp.1d_delt


          python ${QCDIR}/save_image_diff.py -o ${OUTDIR} -i ${OUTDIR} -a r2${PREF}func_data.nii -b rp${PREF}func_data.nii -msg1 'PESTICA+3dvorleg' -msg2 'PESTICA+SLOMOCO'
          cp ${OUTDIR}/slomoco4/*.txt ${OUTDIR}/
          cp ${OUTDIR}/slomoco4/*.1D ${OUTDIR}/
          cp ${OUTDIR}/slomoco4/*.png ${OUTDIR}/

      } &> ${OUTDIR}/logs/moco.log

      if [ ! -f "${OUTDIR}/rp${PREF}func_data.nii" ]; then
          print_time 'ERROR'
          echo "Motion correction did not run successfully."
          echo "Please check ${OUTDIR}/logs/moco.log for more info."
          echo "Execution will abort."
          exit 1
      fi

      PREF=rp${PREF}
fi


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if [ "$DO_SLC" -eq "1" ]; then
      print_time 'INFO'
      echo "Performing slice timing correction using the filtershift method. [${OUTDIR}/logs/slice_timing.log]"
      {
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

          filtershift --in=${OUTDIR}/${PREF}func_data.nii --itl=${ACQ_TYPE}  --cf=${CF} --rt=${REF_TIME} --TR=${TR} --out=${OUTDIR}/a${PREF}func_data.nii
          gunzip -f ${OUTDIR}/a${PREF}func_data.nii.gz

      } &> ${OUTDIR}/logs/slice_timing.log

      if [ ! -f "${OUTDIR}/a${PREF}func_data.nii" ]; then
          print_time 'ERROR'
          echo "Slice timing correction did not run successfully."
          echo "Please check ${OUTDIR}/logs/slice_timing.log for more info."
          echo "Execution will abort."
          exit 1
      fi

      PREF=a${PREF}
fi

if [ "$DO_DEOBL" -eq "1" ]; then
      print_time 'INFO'
      echo "Removing obliqueness information (3dresample + 3drefit). []${OUTDIR}/logs/deoblique.log]"
      {
      3dresample -prefix ${OUTDIR}/w${PREF}func_data.nii -orient lpi -input ${OUTDIR}/${PREF}func_data.nii
      3drefit -deoblique ${OUTDIR}/w${PREF}func_data.nii

    } &> ${OUTDIR}/deoblique.log

    if [ ! -f "${OUTDIR}/w${PREF}func_data.nii" ]; then
        print_time 'ERROR'
        echo "Deoblique did not run successfully."
        echo "Please check ${OUTDIR}/logs/deoblique.log for more info."
        echo "Execution will abort."
        exit 1
    fi

    PREF=w${PREF}
fi


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Skull stripping
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
print_time 'INFO'
echo "Performing brain extraction."
echo "N.B.: This is performed on individual 3d volumes as motion correction might not have been performed."
{
mkdir ${OUTDIR}/tmp
3dTsplit4D -prefix ${OUTDIR}/tmp/func_data.nii -keep_datum ${OUTDIR}/${PREF}func_data.nii
parallel -j4 --line-buffer skull_strip ::: $(seq 0 ${NVOLS}) ::: ${OUTDIR} ::: ${ABIN} ::: ${DO_NLM}

PREF=m${PREF}
3dTcat -prefix ${OUTDIR}/${PREF}func_data.nii -tr ${TR} ${OUTDIR}/tmp/func_data_n*.nii

rm -r ${OUTDIR}/tmp
} &> ${OUTDIR}/logs/masking.log

if [ ! -f "${OUTDIR}/${PREF}func_data.nii" ]; then
    print_time 'ERROR'
    echo "Skull stripping did not run successfully."
    echo "Please check ${OUTDIR}/logs/masking.log for more info."
    echo "Execution will abort."
    exit 1
fi



# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# If movement is extreme, this might help. Then again, it is questionable when exactly to do this.
# At the moment, I'm leaning towards not performing this step and rather censoring offending volumes with the 3dTproject function
if [ "$DO_DPK" -eq "1" ]; then
      print_time 'INFO'
      echo "Removing voxel-wise intensity spikes (3dDespike). [${OUTDIR}/logs/despike.log]"
      {
            3dDespike -nomask -NEW -localedit -prefix ${OUTDIR}/d${PREF}func_data.nii ${OUTDIR}/${PREF}func_data.nii
            python ${QCDIR}/save_image_diff.py -o ${OUTDIR} -i ${OUTDIR} -a ${PREF}func_data.nii -b d${PREF}func_data.nii -msg1 'Before Despike' -msg2 'After Despike'


      } &> ${OUTDIR}/logs/despike.log

      if [ ! -f "${OUTDIR}/d${PREF}func_data.nii" ]; then
          print_time 'ERROR'
          echo "Despike did not run successfully."
          echo "Please check ${OUTDIR}/logs/despike.log for more info."
          echo "Execution will abort."
          exit 1
      fi

      PREF=d${PREF}
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
      print_time 'INFO'
      echo "Performing motion correction using 3dvolreg method [${OUTDIR}/logs/moco.log]"
      {
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

      } &> ${OUTDIR}/logs/moco.log

      if [ ! -f "${OUTDIR}/r${PREF}func_data.nii" ]; then
          print_time 'ERROR'
          echo "Motion correction did not run successfully."
          echo "Please check ${OUTDIR}/logs/moco.log for more info."
          echo "Execution will abort."
          exit 1
      fi

      PREF=r${PREF}
fi


if [ "$DO_DPK" -eq "2" ]; then
      print_time 'INFO'
      echo "Removing voxel-wise intensity spikes (3dDespike). [${OUTDIR}/logs/despike.log]"
      {
          3dDespike -nomask -NEW -localedit -prefix ${OUTDIR}/d${PREF}func_data.nii ${OUTDIR}/${PREF}func_data.nii
          python ${QCDIR}/save_image_diff.py -o ${OUTDIR} -i ${OUTDIR} -a ${PREF}func_data.nii -b d${PREF}func_data.nii -msg1 'Before Despike' -msg2 'After Despike'

      } &> ${OUTDIR}/logs/despike.log

      if [ ! -f "${OUTDIR}/d${PREF}func_data.nii" ]; then
          print_time 'ERROR'
          echo "Despike did not run successfully."
          echo "Please check ${OUTDIR}/logs/despike.log for more info."
          echo "Execution will abort."
          exit 1
      fi

      PREF=d${PREF}
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
      print_time 'INFO'
      echo "Estimating physiological (breathing + heart) noise using PESTICA. [${OUTDIR}/logs/pestica.log]"
      {
      # 0 here indicates that SLOMOCO procedure should not be done
      source ./align_pestica.sh ${OUTDIR} ${PREF}func_data ${TR} 0
      python ${QCDIR}/save_image_diff.py -o ${OUTDIR} -i ${OUTDIR} -a ${PREF}func_data.nii -b p${PREF}func_data.nii -msg1 'Before PESTICA' -msg2 'After PESTICA'

      } &> ${OUTDIR}/logs/pestica.log

      if [ ! -f "${OUTDIR}/p${PREF}func_data.nii" ]; then
          print_time 'ERROR'
          echo "PESTICA did not run successfully."
          echo "Please check ${OUTDIR}/logs/pestica.log for more info."
          echo "Execution will abort."
          exit 1
      fi

      PREF=p${PREF}
fi



# Create native space brain mask
{
3dAutomask -prefix ${OUTDIR}/nat_mask.nii ${OUTDIR}/${PREF}func_data.nii
} &> /dev/null


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Distortion correction using reverse phase acquisition
if [ "$DO_FMAP" -eq "1" ]; then
      print_time 'INFO'
      echo "Geometric distortion correction using TOPUP. [${OUTDIR}/logs/fieldmap.log]"
      {
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


      python ../util/write_topup_file.py -out ${OUTDIR}/topupfield.txt

      # Fieldmap estimation [ONLY estimation]
      topup --imain=${OUTDIR}/topup_data.nii --datain=${OUTDIR}/topupfield.txt --out=${OUTDIR}/topup_results --fout=${OUTDIR}/fieldmap --iout=${OUTDIR}/w_topup_data \
            --estmov=0,0,0  --regmod=membrane_energy  --minmet=1,1,1 --verbose --warpres=8,6,4 --miter=32,12,4 --subsamp=2,2,1 --fwhm=6,4,2

      applytopup --imain=${OUTDIR}/${PREF}func_data.nii --datain=${OUTDIR}/topupfield.txt --topup=${OUTDIR}/topup_results --inindex=1 --method=jac --interp=spline --out=${OUTDIR}/u${PREF}func_data
      gunzip -f ${OUTDIR}/u${PREF}func_data.nii.gz
      3dAutomask -apply_prefix ${OUTDIR}/tmp.nii ${OUTDIR}/u${PREF}func_data.nii
      rm ${OUTDIR}/u${PREF}func_data.nii
      mv ${OUTDIR}/tmp.nii ${OUTDIR}/u${PREF}func_data.nii
      python ${QCDIR}/save_image_diff.py -o ${OUTDIR} -i ${OUTDIR} -a ${PREF}func_data.nii -b u${PREF}func_data.nii -type mean -msg1 'Before FMAP' -msg2 'After FMAP'

      } &> ${OUTDIR}/logs/fieldmap.log


      if [ ! -f "${OUTDIR}/u${PREF}func_data.nii" ]; then
          print_time 'ERROR'
          echo "Geometric distortion correction did not run successfully."
          echo "Please check ${OUTDIR}/logs/fieldmap.log for more info."
          echo "Execution will abort."
          exit 1
      fi
      PREF=u${PREF}
fi


# Distortion correction using acquired fieldmap (magnitude and phase)
if [ "$DO_FMAP" -eq "2" ]; then
      print_time 'INFO'
      echo "Geometric distortion correction using fugue and acquired fieldmap. [${OUTDIR}/logs/fieldmap.log]"
      {
            fslmaths ${OUTDIR}/${PREF}func_data -Tmean  ${OUTDIR}/mean_func


            flirt  -in ${OUTDIR}/fmap_mag_brain -ref ${OUTDIR}/mean_func -omat ${OUTDIR}/fmap2func

            flirt  -in ${OUTDIR}/fmap_mag -ref ${OUTDIR}/mean_func -applyxfm -init ${OUTDIR}/fmap2func -out ${OUTDIR}/fmap_mag_func
            flirt  -in ${OUTDIR}/fmap_phase_rads -ref ${OUTDIR}/mean_func -applyxfm -init ${OUTDIR}/fmap2func -out ${OUTDIR}/fmap_phase_rads_func


            fugue --loadfmap=${OUTDIR}/fmap_phase_rads_func -s 1 --despike --savefmap=${OUTDIR}/fmap_phase_rads_func

            #TODO Pass unwarpdir to configuration
            fugue -i ${OUTDIR}/${PREF}func_data.nii --dwell=${EES} --loadfmap=${OUTDIR}/fmap_phase_rads_func -u ${OUTDIR}/u${PREF}func_data.nii --unwarpdir=y-

            gunzip -f ${OUTDIR}/u${PREF}func_data.nii.gz
            python ${QCDIR}/save_image_diff.py -o ${OUTDIR} -i ${OUTDIR} -a ${PREF}func_data.nii -b u${PREF}func_data.nii -type mean -msg1 'Before FMAP' -msg2 'After FMAP'
      } &> ${OUTDIR}/logs/fieldmap.log


      if [ ! -f "${OUTDIR}/u${PREF}func_data.nii" ]; then
          print_time 'ERROR'
          echo "Geometric distortion correction did not run successfully."
          echo "Please check ${OUTDIR}/logs/fieldmap.log for more info."
          echo "Execution will abort."
          exit 1
      fi
      PREF=u${PREF}
fi


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if [ "$DO_ING" -eq "1" ]; then
      print_time 'INFO'
      echo "Global intensity normalisation to a mean of ${GLOB_VAL}."
      {
        fslmaths ${OUTDIR}/${PREF}func_data -ing ${GLOB_VAL} ${OUTDIR}/i${PREF}func_data
        gunzip -f ${OUTDIR}/i${PREF}func_data.nii.gz
      } &> /dev/null
      PREF=i${PREF}
fi


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Update brain mask
{
rm ${OUTDIR}/nat_mask.nii
3dAutomask -prefix ${OUTDIR}/nat_mask.nii -apply_prefix ${OUTDIR}/proc_data_native.nii ${OUTDIR}/${PREF}func_data.nii
} &> /dev/null

if [ "$DO_CLEAN" -eq "1" ]; then
      print_time 'INFO'
      echo "Removing intermediate files."
      {
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
    } &> /dev/null
fi

print_time 'INFO'
echo "Preprocessing of functional data in native space has been conclude successfully."
exit 0
