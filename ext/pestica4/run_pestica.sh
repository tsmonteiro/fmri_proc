#! /bin/bash

AFNI_BYTEORDER=LSB_FIRST
export AFNI_BYTEORDER

function Usage () {
  cat <<EOF
  PESTICA_v4

  Usage:  run_pestica.sh -d <epi_filename>
 	     -d=dataset: <epi_filename> is the file prefix of
	     the 3D+time dataset to use and correct
               Note: this script will detect suffix
  Opt :  run_pestica.sh -d <epi_filename> -b
             -b=batch mode: turns off interactive popup at stage 3
  Opt :  run_pestica.sh -d <epi_filename> -r <reference_file>
             -r=fMRI mode: load a reference timeseries and include as regressor
	     of no interest to "protect" it from noise removal. See NITRC post
	     on fmri_protection and the background/theoretical considerations
  Opt :  run_pestica.sh -d <epi_filename> -s "1 2"
	     -s=stages: Run stages, from 1-6.  Allows you to run or re-run parts
	     of PESTICA. Give the stages (in increasing order) that you want to
             run. Otherwise, this script will run all stages in order
         stage 1: ICA decomposition
         stage 2: template coregistration and estimation)
         stage 3: cardiac and respiratory responsed vector temporal filtering
         stage 4: RETROICOR with physiological EPI fluctuation
         stage 5: QA

  PMU: Physiologic Monitoring Unit. Use if you have monitored pulse and respiration.
         Currently only Siemens scanner output is supported.
  Opt :  run_pestica.sh -d <epi_filename> -p <Siemens Physio Files Prefix>
             -p=use PMU
	Two types of PMU data format are supported. a) typical Siemens default pmu
        format, e.g. *.ext, *.resp, *.card files. b) CMRR SMS (MB) EPI physio file,
	e.g. *.log
	check rw_pmu_siemens.m file when having a trouble to read pmu.

	Note: If you use the truncated EPI data set as an input, e.g. removal of
              the first or the last 4 volumes, PMU data and input EPI images do not
              have the same length, and you have to inform it correctly here

        a) default;
           rw_pmu_siemens.m assumes the first volume(s) of EPI data is(are) truncated
           as many as the discrepancy of the lengths. Therefore, no input is needed

        b) In case of the last N volume removal of EPI set
          Use -t option to inform how many time points of PMU data will be discard

          run_pestica.sh -d <epi_filename> -p <Siemens Physio Files Prefix> -t N

   NOTE: if you have un-equilibrated volumes at the start, you have to remove them
         before running PESTICA. Most scanners take "dummy" volumes, where the ADCs are
         turned off but the RF and gradients are running as normal, for the first ~3 seconds
         (modulo TR), but in some scanners this is not so and you can see contrast change
         from the 1st to 2nd volumes. Look at your data!

         Since PESTICA runs temporal ICA, un-equilibrated signals generate bias on estimated PMU

         Recommended: 3dvolreg and discard 1st 4:
         3dvolreg -prefix ep2d_pace.moco -base "ep2d_pace[0]" -1Dmatrix_save motion.1D \\
	          -1Dfile motion.txt -zpad 8 -maxite 60 -heptic "ep2d_pace[4..$]"
         This produces a 3D+time dataset that is motion corrected and is 4
         volumes shorter than the original.  PESTICA estimators seem better if
         it is run on the moco'ed data. Next, "run_pestica.sh -d ep2d_pace.moco"

         You can test first volumes for spin saturation with:
           3dToutcount <epi_filename> | 1dplot -stdin -one
         Is the first volume much higher than rest? If so, you may need to remove first
         several volumes first. If you don't know what this means, consult someone who does know,
         this is very important, regression corrections (and analyses) perform poorly
         when the data has unsaturated volumes at the start

EOF
  exit 1
}

allstagesflag=1; batchflag=0; pmuflag=0;  stagepmu1234flag=0;
stage1flag=0; stage2flag=0; stage3flag=0; stage4flag=0; stage5flag=0;
reference="9pestica_null_reference_protect9"
nVolEndCutOff=0  # no EPI volumes at the end were truncated as default

#MATLAB_AFNI_DIR=${MATLAB_AFNI_DIR}

echo nVolEndCutOff=$nVolEndCutOff
while getopts hd:r:fp:bt:s: opt; do
  case $opt in
    h)
       Usage
       exit 1
       ;;
    d) # base 3D+time EPI dataset to use to perform all ICA decomposition and alignment
       epi=$OPTARG
       ;;
    r) # reference timeseries to use as regressor of no interest to "protect" activation
       reference=$OPTARG
       ;;
    s) # option to run stages manually
       allstagesflag=0
       stages=$OPTARG
       # loop over numbers listed
       for i in $stages ; do
         if [ $i -eq 1 ] ; then
           stage1flag=1
	 elif [ $i -eq 2 ] ; then
           stage2flag=1
	 elif [ $i -eq 3 ] ; then
           stage3flag=1
	 elif [ $i -eq 4 ] ; then
           stage4flag=1
	 elif [ $i -eq 5 ] ; then
           stage5flag=1
	 else
	   echo "incorrect syntax for stages input: $i"
	   Usage
	   exit 1
	 fi
       done
      ;;
    b) # flag for batch mode (non-interactive)
       batchflag=1
       ;;
    p) # load monitored pulse ox and respiratory bellows from Siemens PMU system
       pmufileprefix=$OPTARG
       pmuflag=1
       ;;
    t) # the number of volumes truncted at the end of the EPI acquisitions
       nVolEndCutOff=$OPTARG
       echo nVolEndCutOff=$nVolEndCutOff
       ;;
    :)
      echo "option requires input"
      exit 1
      ;;
  esac
done

if [ $pmuflag -eq 1 ]; then
  if [ $allstagesflag -eq 1 ] ; then
    stagepmu1234flag=1; stage5flag=1;
  fi
  if [ $stage1flag -eq 1 ] | [ $stage2flag -eq 1 ] | [ $stage3flag -eq 1 ] | [ $stage4flag -eq 1 ] ; then
    stagepmu1234flag=1;
  fi
else
  if [ $allstagesflag -eq 1 ] ; then
    stage1flag=1; stage2flag=1; stage3flag=1; stage4flag=1; stage5flag=1;
  fi
fi

pesticav=pestica4
if [ $pmuflag -eq 1 ]; then
  pmusuffix="_pmu"
fi

fullcommand="$0"
PESTICA_DIR=`dirname $fullcommand`
export PESTICA_DIR=$PESTICA_DIR

# first test if we are running run_pestica.sh from the base PESTICA_DIR
homedir=`pwd`
if [ $homedir == $PESTICA_DIR ] ; then
  echo "you cannot run PESTICA from the downloaded/extracted PESTICA_DIR"
  echo "please run this from the directory containing the data (or copy of the data)"
  echo "that you want to correct.  Exiting..."
  exit 1
fi

# test for presence of input EPI file with one of the accepted formats and set suffix
if [ -f $epi.hdr ] ; then
  suffix=".hdr"
elif [ -f $epi+orig.HEAD ] ; then
  suffix="+orig.HEAD"
elif [ -f $epi.nii ] ; then
  suffix=".nii"
elif [ -f $epi.nii.gz ] ; then
  suffix=".nii.gz"
else
  echo "3D+time EPI dataset $epi must exist, check filename (do not give suffix)"
  echo "accepted formats: hdr  +orig  nii  nii.gz"
  Usage
  echo ""
  echo "*****   $epi does not exist, exiting ..."
  echo "accepted formats: hdr  +orig  nii  nii.gz"
  exit 1
fi

# create PESTICA subdirectory if it does not exist
epi_pestica="${pesticav}"
if [ ! -d $epi_pestica ] ; then
 echo ""
 echo "* Creating PESTICA Directory: $epi_pestica"
 echo ""
 mkdir $epi_pestica
 echo "mkdir $epi_pestica" >> $epi_pestica/pestica_history.txt
 echo "" >> $epi_pestica/pestica_history.txt
fi

echo "*****   Using $epi+orig.HEAD as input timeseries"
# always copy input file into PESTICA subdirectory in AFNI format
if [ ! -f $epi_pestica/$epi+orig.HEAD ] ; then
  echo "Copying: 3dcopy $epi$suffix $epi_pestica/$epi"
  3dcopy $epi$suffix $epi_pestica/$epi
  echo "3dcopy $epi$suffix $epi_pestica/$epi" >> $epi_pestica/pestica_history.txt
  echo "" >> $epi_pestica/pestica_history.txt
fi

# write command line and PESTICA_DIR to history file
echo "`date`" >> $epi_pestica/pestica_history.txt
echo "`svn info $PESTICA_DIR/run_pestica.sh |grep URL`" >> $epi_pestica/pestica_history.txt
echo "`svn info $PESTICA_DIR/run_pestica.sh |grep Rev`" >> $epi_pestica/pestica_history.txt
echo "PESTICA_v4.0 command line: `basename $fullcommand` $*" >> $epi_pestica/pestica_history.txt
echo "PESTICA env:
`env | grep PESTICA`" >> $epi_pestica/pestica_history.txt
echo "" >> $epi_pestica/pestica_history.txt

# do ALL work inside the pestica/ subdirectory
homedir=`pwd`
cd $epi_pestica
echo "cd $epi_pestica" >> pestica_history.txt
echo "" >> pestica_history.txt
epi_mask="$epi".brain

# test for presence of mask file
if [ ! -f $epi_mask+orig.HEAD ] ; then
  echo ""
  echo "*****   $epi_mask+orig.HEAD does not exist, creating mask"
  echo "note, if you wish to use your own mask/brain file, kill this script"
  echo "then 3dcopy your mask/brain to $epi_mask and re-run."
  echo "PESTICA will use whatever resides in $epi_mask"
  echo ""
  echo "running 3dSkullStrip -input $epi+orig -prefix $epi_mask"
  sleep 1
  3dSkullStrip -input $epi+orig -prefix ___tmp_mask
  #3dAutomask -dilate 4 -prefix $epi_mask $epi+orig
  # erode mask by one voxel
  3dcalc -a ___tmp_mask+orig -prefix ___tmp_mask_ones+orig -expr 'step(a)'
  3dcalc -a ___tmp_mask_ones+orig -prefix ___tmp_mask_ones_dil -b a+i -c a-i -d a+j -e a-j -f a+k -g a-k -expr 'amongst(1,a,b,c,d,e,f,g)'
  3dcalc -a "$epi+orig[0]" -b ___tmp_mask_ones_dil+orig -prefix $epi_mask -expr 'a*step(b)'
  rm ___tmp_mask*
  echo ""
  echo "done with skull-stripping - please check file and if not satisfied, I recommend running"
  echo "3dSkullStrip with different parameters to attempt to get a satisfactory brain mask."
  echo "Either way, this script looks in $epi_pestica/ for $epi_mask to use as your brain mask/strip"
  sleep 1
  echo "3dSkullStrip -input $epi+orig -prefix ___tmp_mask" >> pestica_history.txt
  echo "3dcalc -a ___tmp_mask+orig -prefix ___tmp_mask_ones+orig -expr 'step(a)'" >> pestica_history.txt
  echo "3dcalc -a ___tmp_mask_ones+orig -prefix ___tmp_mask_ones_dil -b a+i -c a-i -d a+j -e a-j -f a+k -g a-k -expr 'amongst(1,a,b,c,d,e,f,g)'" >> pestica_history.txt
  echo "3dcalc -a "$epi+orig[0]" -b ___tmp_mask_ones_dil+orig -prefix $epi_mask -expr 'a*step(b)'" >> pestica_history.txt
  echo "rm ___tmp_mask*" >> pestica_history.txt
  echo "" >> pestica_history.txt
fi
echo "*****   Using $epi_mask+orig.HEAD to mask out non-brain voxels"

### Optional PMU data formatting
if [ $stagepmu1234flag -eq 1 ] ; then
  echo "Using files with prefix: $pmufileprefix"
  pmufileprefix="../$pmufileprefix"
  # convert PhysioLog files into usable data

  echo "matlab $MATLABLINE disp('Starting script...'); addpath $MATLAB_PESTICA_DIR; addpath $MATLAB_AFNI_DIR; rw_pmu_siemens('$epi+orig','$pmufileprefix',$nVolEndCutOff); exit;"
  echo "matlab $MATLABLINE disp('Starting script...'); addpath $MATLAB_PESTICA_DIR; addpath $MATLAB_AFNI_DIR; rw_pmu_siemens('$epi+orig','$pmufileprefix',$nVolEndCutOff); exit;"  >> pestica_history.txt
    matlab $MATLABLINE <<<"disp('Starting script...'); addpath $MATLAB_PESTICA_DIR; addpath $MATLAB_AFNI_DIR; rw_pmu_siemens('$epi+orig','$pmufileprefix',$nVolEndCutOff); exit;"

  # Get IRFs in data - the cardiac estimators are dithered slightly due to the filtering (we will fix that below...)
  echo "Getting IRFs and correcting data with IRF-RETROICOR"
  echo "writing out corrected data as: irfretroicor"
  echo "writing out statistical coupling maps as: coupling_irfret_card coupling_irfret_resp"
  echo "matlab  $MATLABLINE disp('Wait, script starting...'); addpath $MATLAB_PESTICA_DIR; addpath $MATLAB_AFNI_DIR; addpath $MATLAB_EEGLAB_DIR; load RetroTSpmu.mat; [CARD RESP] = retroicor_get_irf('$epi+orig',CARD,RESP,2,'$epi_mask+orig'); irf_retroicor('$epi+orig',CARD,RESP,'$epi_mask+orig','../$reference'); exit;"
  echo "matlab  $MATLABLINE disp('Wait, script starting...'); addpath $MATLAB_PESTICA_DIR; addpath $MATLAB_AFNI_DIR; addpath $MATLAB_EEGLAB_DIR; load RetroTSpmu.mat; [CARD RESP] = retroicor_get_irf('$epi+orig',CARD,RESP,2,'$epi_mask+orig'); irf_retroicor('$epi+orig',CARD,RESP,'$epi_mask+orig','../$reference'); exit;" >> pestica_history.txt
    matlab  $MATLABLINE <<<"disp('Wait, script starting...'); addpath $MATLAB_PESTICA_DIR; addpath $MATLAB_AFNI_DIR; addpath $MATLAB_EEGLAB_DIR; load RetroTSpmu.mat; [CARD RESP] = retroicor_get_irf('$epi+orig',CARD,RESP,2,'$epi_mask+orig'); irf_retroicor('$epi+orig',CARD,RESP,'$epi_mask+orig','../$reference'); exit;"

fi
########### End PMU data correction ###########



########### Start STAGE 1 ###########
if [[ $stage1flag -eq 1 ]] ; then
  echo ""
  echo "Running Stage 1: slicewise temporal Infomax ICA"
  echo ""
        # Use of <<< for MATLAB input was contributed by I. Schwabacher.
  echo "matlab $MATLABLINE addpath $MATLAB_AFNI_DIR; addpath $MATLAB_PESTICA_DIR; addpath $MATLAB_EEGLAB_DIR;disp('Wait, script starting...'); prepare_ICA_decomp(15,'$epi+orig','$epi_mask+orig'); disp('Stage 1 Done!'); exit;"
  echo "matlab $MATLABLINE addpath $MATLAB_AFNI_DIR; addpath $MATLAB_PESTICA_DIR; addpath $MATLAB_EEGLAB_DIR;disp('Wait, script starting...'); prepare_ICA_decomp(15,'$epi+orig','$epi_mask+orig'); disp('Stage 1 Done!'); exit;" >> pestica_history.txt
    matlab $MATLABLINE  <<<"addpath $MATLAB_AFNI_DIR; addpath $MATLAB_PESTICA_DIR; addpath $MATLAB_EEGLAB_DIR;disp('Wait, script starting...'); prepare_ICA_decomp(15,'$epi+orig','$epi_mask+orig'); disp('Stage 1 Done!'); exit;"
fi
########### End STAGE 1 ###########

########### Start STAGE 2 ###########
if [[ $stage2flag -eq 1 ]] ; then
  echo ""
  echo "Running Stage 2: Coregistration of EPI to MNI space and back-transform of templates, followed by PESTICA estimation"
  echo ""
  ## EPI to MNI
  #3dcopy $PESTICA_VOL_DIR/resp_mean_mni_${pesticav}.nii $PESTICA_VOL_DIR/resp_mean_mni_${pesticav}+orig
  #3dcopy $PESTICA_VOL_DIR/card_mean_mni_${pesticav}.nii $PESTICA_VOL_DIR/card_mean_mni_${pesticav}+orig
  #3dcopy $PESTICA_VOL_DIR/meanepi_mni.nii $PESTICA_VOL_DIR/meanepi_mni+orig

  #3drefit -byteorder MSB_FIRST  $PESTICA_VOL_DIR/resp_mean_mni_${pesticav}+orig
  #3drefit -byteorder MSB_FIRST  $PESTICA_VOL_DIR/card_mean_mni_${pesticav}+orig
  #3drefit -byteorder MSB_FIRST  $PESTICA_VOL_DIR/meanepi_mni+orig

  if [ ! -f mni.coreg.$epi_mask.1D ] ; then

    echo "Coregistration to EPI template"
    echo 3dAllineate -automask -prefix ./$epi_mask.crg2mni.nii -source $epi_mask+orig -base $PESTICA_VOL_DIR/meanepi_mni.nii -1Dmatrix_save $epi_mask.coreg.mni.1D
    echo "3dAllineate -automask -prefix ./$epi_mask.crg2mni.nii -source $epi_mask+orig -base $PESTICA_VOL_DIR/meanepi_mni.nii -1Dmatrix_save $epi_mask.coreg.mni.1D" >> pestica_history.txt
          3dAllineate -nomask  -prefix ./$epi_mask.crg2mni.nii -source $epi_mask+orig -base $PESTICA_VOL_DIR/meanepi_mni.nii -1Dmatrix_save $epi_mask.coreg.mni.1D

         cat_matvec $epi_mask.coreg.mni.1D -I -ONELINE > mni.coreg.$epi_mask.1D

#    3dcopy $epi_mask.crg2mni.nii $epi_mask.crg2mni+orig
#    3drefit -byteorder MSB_FIRST  $epi_mask.crg2mni+orig

  fi



  # move PESTICA template mni to EPI space
  echo "3dAllineate -prefix ./resp_${pesticav}.nii -source $PESTICA_VOL_DIR/resp_mean_mni_${pesticav}.nii -base $epi_mask+orig -1Dmatrix_apply mni.coreg.$epi_mask.1D -overwrite"
  echo "3dAllineate -prefix ./resp_${pesticav}.nii -source $PESTICA_VOL_DIR/resp_mean_mni_${pesticav}.nii -base $epi_mask+orig -1Dmatrix_apply mni.coreg.$epi_mask.1D -overwrite" >> pestica_history.txt
        3dAllineate -prefix ./resp_${pesticav}.nii -source $PESTICA_VOL_DIR/resp_mean_mni_${pesticav}.nii -base $epi_mask+orig -1Dmatrix_apply mni.coreg.$epi_mask.1D -overwrite
  echo "3dAllineate -prefix ./card_${pesticav}.nii -source $PESTICA_VOL_DIR/card_mean_mni_${pesticav}.nii -base $epi_mask+orig -1Dmatrix_apply mni.coreg.$epi_mask.1D -overwrite"
  echo "3dAllineate -prefix ./card_${pesticav}.nii -source $PESTICA_VOL_DIR/card_mean_mni_${pesticav}.nii -base $epi_mask+orig -1Dmatrix_apply mni.coreg.$epi_mask.1D -overwrite" >> pestica_history.txt
        3dAllineate -prefix ./card_${pesticav}.nii -source $PESTICA_VOL_DIR/card_mean_mni_${pesticav}.nii -base $epi_mask+orig -1Dmatrix_apply mni.coreg.$epi_mask.1D -overwrite

# 3dcopy card_${pesticav}.nii card_${pesticav}+orig
# 3dcopy resp_${pesticav}.nii resp_${pesticav}+orig


# 3drefit -byteorder MSB_FIRST  card_${pesticav}+orig
# 3drefit -byteorder MSB_FIRST  resp_${pesticav}+orig

  # run PESTICA
  echo "Obtaining PESTICA estimators"
  echo  "matlab $MATLABLINE addpath $MATLAB_AFNI_DIR; addpath $MATLAB_PESTICA_DIR; [card,resp]=apply_PESTICA(15,'$epi+orig','$epi_mask+orig','${pesticav}'); fp=fopen('card_raw_${pesticav}.dat','w'); fprintf(fp,'%g\n',card); fclose(fp); fp=fopen('resp_raw_${pesticav}.dat','w'); fprintf(fp,'%g\n',resp); fclose(fp); exit"
  echo  "matlab $MATLABLINE addpath $MATLAB_AFNI_DIR; addpath $MATLAB_PESTICA_DIR; [card,resp]=apply_PESTICA(15,'$epi+orig','$epi_mask+orig','${pesticav}'); fp=fopen('card_raw_${pesticav}.dat','w'); fprintf(fp,'%g\n',card); fclose(fp); fp=fopen('resp_raw_${pesticav}.dat','w'); fprintf(fp,'%g\n',resp); fclose(fp); exit" >> pestica_history.txt
         matlab $MATLABLINE <<<"addpath $MATLAB_AFNI_DIR; addpath $MATLAB_PESTICA_DIR; disp('Wait, script starting...'); [card,resp]=apply_PESTICA(15,'$epi+orig','$epi_mask+orig','${pesticav}'); fp=fopen('card_raw_${pesticav}.dat','w'); fprintf(fp,'%g\n',card); fclose(fp); fp=fopen('resp_raw_${pesticav}.dat','w'); fprintf(fp,'%g\n',resp); fclose(fp); disp('Stage 2 Done!'); exit;"
fi
########### End STAGE 2 ###########

########### Start STAGE 3 ###########
if [[ $stage3flag -eq 1 ]] ; then
  echo ""
  echo "Running Stage 3: Filtering PESTICA estimators, cardiac first, then respiratory"
  echo ""
  echo "NOTE: TR must be set correctly in header for 3D+time dataset - if in doubt, check it and correct it"
  echo "Values given below:"
  echo "matlab $MATLABLINE addpath $MATLAB_PESTICA_DIR; addpath $MATLAB_AFNI_DIR; load('card_raw_${pesticav}.dat'); load('resp_raw_${pesticav}.dat'); card=view_and_correct_estimator(card_raw_${pesticav},'$epi+orig','c',$batchflag); resp=view_and_correct_estimator(resp_raw_${pesticav},'$epi+orig','r',$batchflag); fp=fopen('card_${pesticav}.dat','w'); fprintf(fp,'%g\n',card); fclose(fp); fp=fopen('resp_${pesticav}.dat','w'); fprintf(fp,'%g\n',resp); fclose(fp); exit;"
  echo "matlab $MATLABLINE addpath $MATLAB_PESTICA_DIR; addpath $MATLAB_AFNI_DIR;load('card_raw_${pesticav}.dat'); load('resp_raw_${pesticav}.dat'); card=view_and_correct_estimator(card_raw_${pesticav},'$epi+orig','c',$batchflag); resp=view_and_correct_estimator(resp_raw_${pesticav},'$epi+orig','r',$batchflag); fp=fopen('card_${pesticav}.dat','w'); fprintf(fp,'%g\n',card); fclose(fp); fp=fopen('resp_${pesticav}.dat','w'); fprintf(fp,'%g\n',resp); fclose(fp); exit;" >> pestica_history.txt
        matlab $MATLABLINE <<<"addpath $MATLAB_PESTICA_DIR; addpath $MATLAB_AFNI_DIR;load('card_raw_${pesticav}.dat'); load('resp_raw_${pesticav}.dat'); disp('Wait, script starting...'); card=view_and_correct_estimator(card_raw_${pesticav},'$epi+orig','c',$batchflag); resp=view_and_correct_estimator(resp_raw_${pesticav},'$epi+orig','r',$batchflag);  fp=fopen('card_${pesticav}.dat','w'); fprintf(fp,'%g\n',card); fclose(fp); fp=fopen('resp_${pesticav}.dat','w'); fprintf(fp,'%g\n',resp); fclose(fp); disp('Stage 3 Done!'); exit;"
fi
########### End STAGE 3 ##########

########### Start STAGE 4 ###########
if [[ $stage4flag -eq 1 ]] ; then
# PESTICA convert to phase using 3dretroicor
  echo "Running MATLAB-version of RETROICOR with physiological noise fluctuation"
  echo "matlab  $MATLABLINE addpath $PESTICA_RETROICOR_DIR; addpath $MATLAB_AFNI_DIR;card=textread('card_${pesticav}.dat'); resp=textread('resp_${pesticav}.dat'); retroicor_pestica('$epi+orig',card,resp,2,'$epi_mask+orig'); exit"
  echo "matlab  $MATLABLINE addpath $PESTICA_RETROICOR_DIR; addpath $MATLAB_AFNI_DIR;card=textread('card_${pesticav}.dat'); resp=textread('resp_${pesticav}.dat'); retroicor_pestica('$epi+orig',card,resp,2,'$epi_mask+orig'); exit" >> pestica_history.txt
          matlab  $MATLABLINE <<<"addpath $MATLAB_PESTICA_DIR; addpath $MATLAB_AFNI_DIR;card=textread('card_${pesticav}.dat'); resp=textread('resp_${pesticav}.dat'); disp('Wait, script starting...'); retroicor_pestica('$epi+orig',card,resp,2,'$epi_mask+orig'); disp('Stage 4 done!'); exit;"
fi
########### End STAGE 4 ###########

########### Start STAGE 5 ###########
if [[ $stage5flag -eq 1 ]] ; then
  echo ""
  echo "Running Stage 5: Make QA plots"
  echo ""

  # run QA on output files
  echo "matlab $MATLABLINE addpath $MATLAB_PESTICA_DIR; physio_qa('$epi+orig',$pmuflag);"
  echo "matlab $MATLABLINE addpath $MATLAB_PESTICA_DIR; physio_qa('$epi+orig',$pmuflag);" >> pestica_history.txt
    matlab $MATLABLINE <<<"addpath $MATLAB_PESTICA_DIR; addpath $MATLAB_AFNI_DIR; disp('Wait, script starting...'); physio_qa('$epi+orig',$pmuflag); disp('Done!'); exit;"

  echo " *******************************************************"
  echo " *******************************************************"
  echo " WARNING, AFNI IS ABOUT TO STEAL WINDOW FOCUS!!"
  echo " wait til this script ends in a few seconds, it will end at same time as last AFNI ends"
  echo " sorry, I do not yet have a workaround for AFNI's always stealing focus behavior"
  echo " please email me if you have one"
  echo " *******************************************************"
  echo " *******************************************************"
  sleep 1

  # change this if the plots always give a poor view of the slices - slice 20 in AFNI is reasonable for most acquisitions
  dims=(`3dAttribute DATASET_DIMENSIONS $epi+orig`)
  zdim=${dims[2]}
  if [ $zdim -gt 70 ]; then
    zpos=40
    montstr='6x6'
  elif [ $zdim -gt 50 ]; then
    zpos=30
    montstr='5x5'
  else
    zpos=20
    montstr='4x4'
  fi

  fname=`basename $epi`
  # threshold for cardiac/respiratory coupling is ideally detected from the data itself, but may have to be adjusted manually
  if [ $pmuflag -eq 1 ] ; then
    inamec=coupling_irfret_card_pmu
    inamer=coupling_irfret_resp_pmu
    snamec=pmu_card_coupling_overlay
    snamer=pmu_resp_coupling_overlay
  else
    inamec=coupling_ret_card_pestica
    inamer=coupling_ret_resp_pestica
    snamec=pestica_card_coupling_overlay
    snamer=pestica_resp_coupling_overlay
  fi

  if [ $pmuflag -eq 1 ] ; then
    thresh=`head -n 1 coupling_pmu_thresholds.txt`
    threshr=`tail -n 1 coupling_pmu_thresholds.txt`
  else
    thresh=`head -n 1 coupling_pestica_thresholds.txt`
    threshr=`tail -n 1 coupling_pestica_thresholds.txt`
  fi

    afni -com "OPEN_WINDOW A.axialimage mont="$montstr":2:0:none opacity=6" \
         -com "SET_THRESHOLD A.$thresh 1"  -com 'SET_PBAR_NUMBER A.12' \
         -com 'SET_FUNC_RESAM A.Cu.Cu' -com "SET_UNDERLAY A.$fname" -com 'SET_FUNC_RANGE A.10' \
         -com "SET_OVERLAY A."$inamec" 1 0"        -com 'OPEN_PANEL A.Define_Overlay' -com 'SET_FUNC_VISIBLE A.+' \
         -com "SET_IJK A 64 64 $zpos" -com 'SET_XHAIRS A.OFF' -com "SAVE_PNG A.axialimage $snamec" \
         -com "OPEN_WINDOW A.axialimage mont="$montstr":2:0:none opacity=6" \
         -com "SET_THRESHOLD A.$threshr 1"  -com 'SET_PBAR_NUMBER A.12' \
         -com 'SET_FUNC_RESAM A.Cu.Cu' -com "SET_UNDERLAY A.$fname" -com 'SET_FUNC_RANGE A.10' \
         -com "SET_OVERLAY A."$inamer" 1 0"        -com 'OPEN_PANEL A.Define_Overlay' -com 'SET_FUNC_VISIBLE A.+' \
         -com "SET_IJK A 64 64 $zpos" -com 'SET_XHAIRS A.OFF' -com "SAVE_PNG A.axialimage $snamer" \
         -com 'QUIT' -dset $fname+orig $inamec+orig $inamer+orig >> afnilogfile.txt 2>&1

  # remove copied EPI file inside PESTICA subdir, as we should be finished and don't need to take up the extra space
  rm $epi+orig.????
  echo "rm $epi+orig.???? (temp file removal inside $epi_pestica only)" >> pestica_history.txt
  echo "" >> pestica_history.txt
fi

cd $homedir
echo "End of PESTICA script" >> $epi_pestica/pestica_history.txt
echo "`date`" >> $epi_pestica/pestica_history.txt
echo "" >> $epi_pestica/pestica_history.txt

# Note from Isaac:
# This is necessary because MATLAB likes to break the terminal for some reason.
# but it's commented here because it doesn't like to run in the background.
#reset
