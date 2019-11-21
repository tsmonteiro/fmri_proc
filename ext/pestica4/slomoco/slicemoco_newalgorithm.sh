#! /bin/bash

AFNI_BYTEORDER=LSB_FIRST	
export AFNI_BYTEORDER

function Usage () {
  cat <<EOF
   slicemoco_newalgorithm (distributed with PESTICA v4)
     Algorithm: this script first runs slicewise in-plane (xy) 3DOF motion correction
	        then runs a slicewise 6DOF rigid-body correction for each slice
		this script reads these motion parameters back in and regress on voxel timeseries
       WARNING, make sure you have removed unsaturated images at start (such as first 4 vols)
       You can test first volumes for spin saturation with: 3dToutcount <epi_filename> | 1dplot -stdin -one
       Is the first volume much higher than rest? If so, you may need to remove first several volumes first
       If you don't know what this means, consult someone who does know, this is very important,
       regression corrections (and analyses) perform poorly when the data has unsaturated volumes at the start
        Simple method to remove 1st 4: 3dcalc -a "<epi_file>+orig[4..$]" -expr a -prefix <epi_file>.steadystate

 Usage:  slicemoco_newalgorithm.sh -d <epi_filename>  
 	     -d = dataset: <epi_filename> is the file prefix of
	     the 3D+time dataset to use and correct.
               Note: this script will detect suffix for epi_filename

         slicemoco_newalgorithm.sh -d <epi_filename> -r
             -r = perform in parallel with final PESTICA regression correction
                  this assumes PESTICA estimation steps 1-5 have been run and exist in subdir pestica4/

         slicemoco_newalgorithm.sh -d <epi_filename> -p
             -p = same as -r option, but assuming you used PMU data instead of PESTICA for the correction

 Recommended, run after running PESTICA or PMU correction, so we can incorporate all regressions in parallel:
	       slicemoco_newalgorithm.sh -d <epi_file> -r
	   OR, slicemoco_newalgorithm.sh -d <epi_file> -p

EOF
  exit 1
}

phypest=0
phypmu=0
pesticav=pestica4
slomocov=slomoco4

while getopts hd:pr opt; do
  case $opt in
    h)
       Usage
       exit 1
       ;;
    d) # base 3D+time EPI dataset to use to perform corrections
       epi=$OPTARG
       ;;
    r)
       phypest=1
       ;;
    p)
       phypmu=1
       ;;
    :)
      echo "option requires input"
      exit 1
      ;;
  esac
done

if [ $phypest -eq 1 ] && [ $phypmu -eq 1 ]; then
  echo "cannot set physio correction to both PMU and PESTICA, choose the appropriate one"
  exit 1
elif [ $phypest -eq 0 ] && [ $phypmu -eq 0 ]; then
  echo "Note that the physiologic fluctuation will not be regressed out (not recommended)"
elif [ $phypest -eq 1 ]; then
  echo "Estimated PMU signal will be regressed out."
else
  echo "Measured PMU signal will be regressed out."
fi

# test if epi filename was set
if [ -z $epi ] ; then
  echo "3D+time EPI dataset $epi must be given"
  Usage
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

fullcommand="$0"
if [ -z $SLCMOCO_DIR ] ; then
  echo "setting SLCMOCO_DIR to directory where this script is located (not the CWD)"
  SLCMOCO_DIR=`dirname $fullcommand`
  # this command resides in the moco/ subdirectory of PESTICA_DIR
  export SLCMOCO_DIR=$SLCMOCO_DIR
fi

# first test if we are running run_pestica.sh from the base SLCMOCO_DIR
homedir=`pwd`
if [ $homedir == $SLCMOCO_DIR ] ; then
  echo "you cannot run PESTICA from the downloaded/extracted SLCMOCO_DIR"
  echo "please run this from the directory containing the data (or copy of the data)"
  echo "that you want to correct.  Exiting..."
  exit 1
fi

# create PESTICA subdirectory if it does not exist
epi_pestica=${pesticav}
epi_slomoco=${slomocov}
if [ ! -d $epi_slomoco ] ; then
 echo ""
 echo "* Creating SLOMOCO Directory: $epi_slomoco"
 echo ""
 mkdir $epi_slomoco
 echo "mkdir $epi_slomoco" >> $epi_slomoco/slomoco_history.txt
 echo "" >> $epi_slomoco/slomoco_history.txt
fi

echo "*****   Using copied input in $epi_slomoco/$epi+orig.HEAD as input timeseries"
# always copy input file into SLOMOCO subdirectory in AFNI format
if [ ! -f $epi_slomoco/$epi+orig.HEAD ] ; then
  echo "Copying: 3dcopy $epi$suffix $epi_slomoco/$epi"
  3dcopy $epi$suffix $epi_slomoco/$epi
  echo "3dcopy $epi$suffix $epi_slomoco/$epi" >> $epi_slomoco/slomoco_history.txt
  echo "" >> $epi_slomoco/slomoco_history.txt
fi

# write command line and SLCMOCO_DIR to history file
echo "`date`" >> $epi_slomoco/slomoco_history.txt
echo "PESTICA_afni-v${pesticav} command line: `basename $fullcommand` $*" >> $epi_slomoco/slomoco_history.txt
echo "PESTICA env:
`env | grep PESTICA`" >> $epi_slomoco/slomoco_history.txt
echo "" >> $epi_slomoco/slomoco_history.txt

cd $epi_slomoco
echo "cd $epi_slomoco" >> slomoco_history.txt
echo "" >> slomoco_history.txt

# do ALL work inside the slomoco/ subdirectory
homedir=`pwd`
epi_mask="$epi".brain
if [ ! -f $epi_mask+orig.HEAD ] ; then
  if [ -f ../$epi_pestica/$epi_mask+orig.HEAD ] ; then
    3dcopy ../$epi_pestica/$epi_mask+orig $epi_mask 
  else
    echo ""
    echo "*****   $epi_mask+orig.HEAD does not exist, creating mask"
    echo "note, if you wish to use your own mask/brain file, kill this script"
    echo "then 3dcopy your mask/brain to $epi_mask and re-run."
    echo "SLOMOCO will use whatever resides in $epi_mask"
    echo ""
    echo "running 3dSkullStrip -input $epi+orig -prefix $epi_mask"
    sleep 1
    3dSkullStrip -input $epi+orig -prefix ___tmp_mask
    # erode mask by one voxel
    3dcalc -a ___tmp_mask+orig -prefix ___tmp_mask_ones+orig -expr 'step(a)'
    3dcalc -a ___tmp_mask_ones+orig -prefix ___tmp_mask_ones_dil -b a+i -c a-i -d a+j -e a-j -f a+k -g a-k -expr 'amongst(1,a,b,c,d,e,f,g)'
    3dcalc -a "$epi+orig[0]" -b ___tmp_mask_ones_dil+orig -prefix $epi_mask -expr 'a*step(b)'
    rm ___tmp_mask*
    echo ""
    echo "done with skull-stripping - please check file and if not satisfied, I recommend running"
    echo "3dSkullStrip with different parameters to attempt to get a satisfactory brain mask."
    echo "Either way, this script looks in $epi_slomoco/ for $epi_mask to use as your brain mask/strip"
    echo "      if 3dSkullStrip dies due to large size or other unknown reason, try: "
    echo "   3dAutomask -dilate 1 -clfrac 0.4 -prefix $epi_slomoco/$epi.brain $epi+orig"
    sleep 1
    echo "3dSkullStrip -input $epi+orig -prefix ___tmp_mask" >> slomoco_history.txt
    echo "3dcalc -a ___tmp_mask+orig -prefix ___tmp_mask_ones+orig -expr 'step(a)'" >> slomoco_history.txt
    echo "3dcalc -a ___tmp_mask_ones+orig -prefix ___tmp_mask_ones_dil -b a+i -c a-i -d a+j -e a-j -f a+k -g a-k -expr 'amongst(1,a,b,c,d,e,f,g)'" >> slomoco_history.txt
    echo "3dcalc -a "$epi+orig[0]" -b ___tmp_mask_ones_dil+orig -prefix $epi_mask -expr 'a*step(b)'" >> slomoco_history.txt
    echo "rm ___tmp_mask*" >> slomoco_history.txt
    echo "" >> slomoco_history.txt
  fi
fi
echo "*****   Using $epi_mask+orig.HEAD to mask out non-brain voxels"

# will work regardless of conversion to milliseconds or seconds
SLOMOCO_SLICE_TIMING="`3dAttribute TAXIS_OFFSETS $epi+orig`"
echo $SLOMOCO_SLICE_TIMING > slice_timing.txt    # second unit
SMSfactor=0
if [ "`3dAttribute TAXIS_NUMS $epi+orig | awk '{print $3}'`" -eq 77002 ] ; then
  a=""
  for f in $SLOMOCO_SLICE_TIMING ; do
    a="$a `echo "($f * 1000)" | bc`"
    if [ $f == 0 ]; then
      let "SMSfactor+=1"
    fi
  done
  SLOMOCO_SLICE_TIMING=$a
fi
echo "SLOMOCO_SLICE_TIMING = $SLOMOCO_SLICE_TIMING " >> slomoco_history.txt ## milisecond unit

## SMS acquisition check here ##
if [ $SMSfactor -gt 1 ]; then
  echo SMS acquisition is assumed with MB factor = $SMSfactor
fi

##### calculate vol moco parameter###
if [ -f $epi.mocoafni.txt ] && [ -f $epi.mocoafni.1D ] && [ -f $epi.mocoafni.maxdisp.1D ] ; then
  echo "SKIP; 3dvolume registration has been done."
else
  if [ -f ../$epi.mocoafni.txt ] && [ -f ../$epi.mocoafni.1D ] && [ -f ../$epi.mocoafni.maxdisp.1D ] ; then
    echo ln -s ../$epi.mocoafni.txt $epi.mocoafni.txt
         ln -s ../$epi.mocoafni.txt $epi.mocoafni.txt
    echo ln -s ../$epi.mocoafni.1D  $epi.mocoafni.1D
         ln -s ../$epi.mocoafni.1D  $epi.mocoafni.1D
    echo ln -s ../$epi.mocoafni.maxdisp.1D   $epi.mocoafni.maxdisp.1D 
         ln -s ../$epi.mocoafni.maxdisp.1D   $epi.mocoafni.maxdisp.1D 
  else
    echo 3dvolreg -prefix $epi.mocoafni.hdr -base 0 -zpad 8 -maxite 60 -x_thresh 0.005 -rot_thresh 0.008 -verbose -dfile $epi.mocoafni.txt -1Dfile $epi.mocoafni.1D -maxdisp1D $epi.mocoafni.maxdisp.1D -heptic $epi+orig
    echo "3dvolreg -prefix $epi.mocoafni.hdr -base 0 -zpad 8 -maxite 60 -x_thresh 0.005 -rot_thresh 0.008 -verbose -dfile $epi.mocoafni.txt -1Dfile $epi.mocoafni.1D -maxdisp1D $epi.mocoafni.maxdisp.1D -heptic $epi+orig" >> slomoco_history.txt
         3dvolreg -prefix $epi.mocoafni.hdr -base 0 -zpad 8 -maxite 60 -x_thresh 0.005 -rot_thresh 0.008 -verbose -dfile $epi.mocoafni.txt -1Dfile $epi.mocoafni.1D -maxdisp1D $epi.mocoafni.maxdisp.1D -heptic $epi+orig
    echo "rm $epi.mocoafni.hdr $epi.mocoafni.img" >> slomoco_history.txt
          rm $epi.mocoafni.hdr $epi.mocoafni.img
  fi
  echo "cp $epi.mocoafni.txt mocoafni.txt" >> slomoco_history.txt
        cp $epi.mocoafni.txt mocoafni.txt
fi


if [ ! -d tempslmocoxy_afni_$epi ] ; then
  echo "$SLCMOCO_DIR/run_correction_slicemocoxy_afni.sh -b $epi -p $epi.slicemocoxy_afni"
  echo "$SLCMOCO_DIR/run_correction_slicemocoxy_afni.sh -b $epi -p $epi.slicemocoxy_afni" >> slomoco_history.txt
        $SLCMOCO_DIR/run_correction_slicemocoxy_afni.sh -b $epi -p $epi.slicemocoxy_afni
fi

if [ ! -d tempslmoco_volslc_alg_vol_$epi.slicemocoxy_afni ] ; then
  echo "$SLCMOCO_DIR/run_slicemoco_inside_fixed_vol.sh -b $epi.slicemocoxy_afni" 
  echo "$SLCMOCO_DIR/run_slicemoco_inside_fixed_vol.sh -b $epi.slicemocoxy_afni" >> slomoco_history.txt
        $SLCMOCO_DIR/run_slicemoco_inside_fixed_vol.sh -b $epi.slicemocoxy_afni
fi

##### done vol moco parameter###
echo ""
echo "Running Secondorder Motion Correction using SLOMOCO output"
echo ""
if [ $phypmu -eq 1 ] ; then
  if [ -f ../$epi_pestica/impulse_responses.mat ]; then
    echo "cp ../$epi_pestica/impulse_responses.mat impulse_responses.mat"
    echo "cp ../$epi_pestica/impulse_responses.mat impulse_responses.mat" >> slomoco_history.txt
          cp ../$epi_pestica/impulse_responses.mat impulse_responses.mat

    echo "matlab $MATLABLINE addpath $SLCMOCO_DIR; addpath $MATLAB_AFNI_DIR; addpath $MATLAB_PESTICA_DIR; load impulse_responses.mat; slicemoco_newalgorithm_input('$epi.slicemocoxy_afni+orig','$epi_mask+orig','tempslmoco_volslc_alg_vol_$epi.slicemocoxy_afni/motion.wholevol_zt',CARD,RESP); exit;"
    echo "matlab $MATLABLINE addpath $SLCMOCO_DIR; addpath $MATLAB_AFNI_DIR; addpath $MATLAB_PESTICA_DIR; load impulse_responses.mat; slicemoco_newalgorithm_input('$epi.slicemocoxy_afni+orig','$epi_mask+orig','tempslmoco_volslc_alg_vol_$epi.slicemocoxy_afni/motion.wholevol_zt',CARD,RESP); exit;" >> slomoco_history.txt
          matlab  $MATLABLINE <<<"addpath $SLCMOCO_DIR; addpath $MATLAB_AFNI_DIR; addpath $MATLAB_PESTICA_DIR; load impulse_responses.mat; slicemoco_newalgorithm_input('$epi.slicemocoxy_afni+orig','$epi_mask+orig','tempslmoco_volslc_alg_vol_$epi.slicemocoxy_afni/motion.wholevol_zt',CARD,RESP); exit;"
  corrstr="slicemocoxy_afni.slomoco_pmu"
  else
    echo "ERROR: impulse response mat file is needed"
    exit
  fi
elif [ $phypest -eq 1 ] ; then
  if [ -f ../$epi_pestica/card_${pesticav}.dat ] && [ -f ../$epi_pestica/resp_${pesticav}.dat ] ; then
    echo "cp ../$epi_pestica/card_${pesticav}.dat card_${pesticav}.dat"
    echo "cp ../$epi_pestica/card_${pesticav}.dat card_${pesticav}.dat" >> slomoco_history.txt
          cp ../$epi_pestica/card_${pesticav}.dat card_${pesticav}.dat
    echo "cp ../$epi_pestica/resp_${pesticav}.dat resp_${pesticav}.dat"
    echo "cp ../$epi_pestica/resp_${pesticav}.dat resp_${pesticav}.dat" >> slomoco_history.txt
          cp ../$epi_pestica/resp_${pesticav}.dat resp_${pesticav}.dat

    echo "matlab $MATLABLINE addpath $SLCMOCO_DIR; addpath $MATLAB_AFNI_DIR; addpath $MATLAB_PESTICA_DIR; CARD=load('card_${pesticav}.dat'); RESP=load('resp_${pesticav}.dat'); slicemoco_newalgorithm_input('$epi.slicemocoxy_afni+orig','$epi_mask+orig','tempslmoco_volslc_alg_vol_$epi.slicemocoxy_afni/motion.wholevol_zt',CARD,RESP); exit;"
    echo "matlab $MATLABLINE addpath $SLCMOCO_DIR; addpath $MATLAB_AFNI_DIR; addpath $MATLAB_PESTICA_DIR; CARD=load('card_${pesticav}.dat'); RESP=load('resp_${pesticav}.dat'); slicemoco_newalgorithm_input('$epi.slicemocoxy_afni+orig','$epi_mask+orig','tempslmoco_volslc_alg_vol_$epi.slicemocoxy_afni/motion.wholevol_zt',CARD,RESP); exit;" >> slomoco_history.txt
      matlab $MATLABLINE <<<"addpath $SLCMOCO_DIR; addpath $MATLAB_AFNI_DIR; addpath $MATLAB_PESTICA_DIR; CARD=load('card_${pesticav}.dat'); RESP=load('resp_${pesticav}.dat'); slicemoco_newalgorithm_input('$epi.slicemocoxy_afni+orig','$epi_mask+orig','tempslmoco_volslc_alg_vol_$epi.slicemocoxy_afni/motion.wholevol_zt',CARD,RESP); exit;"
  corrstr="slicemocoxy_afni.slomoco_pestica"
  else
    echo "Error: card_${pesticav}.dat and resp_${pesticav}.dat files should be given"
    exit
  fi
else
  echo "matlab $MATLABLINE addpath $SLCMOCO_DIR; addpath $MATLAB_AFNI_DIR; addpath $MATLAB_PESTICA_DIR; slicemoco_newalgorithm_input('$epi.slicemocoxy_afni+orig','$epi_mask+orig','tempslmoco_volslc_alg_vol_$epi.slicemocoxy_afni/motion.wholevol_zt'); exit;"
  echo "matlab $MATLABLINE addpath $SLCMOCO_DIR; addpath $MATLAB_AFNI_DIR; addpath $MATLAB_PESTICA_DIR; slicemoco_newalgorithm_input('$epi.slicemocoxy_afni+orig','$epi_mask+orig','tempslmoco_volslc_alg_vol_$epi.slicemocoxy_afni/motion.wholevol_zt'); exit;" >> slomoco_history.txt
    matlab $MATLABLINE <<<"addpath $SLCMOCO_DIR; addpath $MATLAB_AFNI_DIR; addpath $MATLAB_PESTICA_DIR; slicemoco_newalgorithm_input('$epi.slicemocoxy_afni+orig','$epi_mask+orig','tempslmoco_volslc_alg_vol_$epi.slicemocoxy_afni/motion.wholevol_zt'); exit;"
  corrstr="slicemocoxy_afni.slomoco"
fi

echo "matlab $MATLABLINE addpath $SLCMOCO_DIR; addpath $MATLAB_PESTICA_DIR; addpath $MATLAB_AFNI_DIR; qa_slomoco('$epi.slicemocoxy_afni+orig','tempslmoco_volslc_alg_vol_$epi.slicemocoxy_afni/motion.wholevol_zt','tempslmocoxy_afni_$epi'); exit;"
echo "matlab $MATLABLINE addpath $SLCMOCO_DIR; addpath $MATLAB_PESTICA_DIR; addpath $MATLAB_AFNI_DIR; qa_slomoco('$epi.slicemocoxy_afni+orig','tempslmoco_volslc_alg_vol_$epi.slicemocoxy_afni/motion.wholevol_zt','tempslmocoxy_afni_$epi'); exit;" >> slomoco_history.txt
  matlab $MATLABLINE <<<"addpath $SLCMOCO_DIR; addpath $MATLAB_PESTICA_DIR; addpath $MATLAB_AFNI_DIR; qa_slomoco('$epi.slicemocoxy_afni+orig','tempslmoco_volslc_alg_vol_$epi.slicemocoxy_afni/motion.wholevol_zt','tempslmocoxy_afni_$epi'); exit;"

# change attributes to match original data
taxis_nums=`3dAttribute TAXIS_NUMS $epi+orig`
taxis_floats=`3dAttribute TAXIS_FLOATS $epi+orig`
taxis_offset=`3dAttribute TAXIS_OFFSETS $epi+orig`

#change TR back to input dataset TR
#3drefit -TR $tr $epi+orig $epi.slicemocoxy_afni.zalg_moco2+orig
#copy over and save t-axis nums into new image registered dataset
echo 3drefit -saveatr -atrint TAXIS_NUMS "$taxis_nums" $epi.$corrstr+orig
     3drefit -saveatr -atrint TAXIS_NUMS "$taxis_nums" $epi.$corrstr+orig
#copy over and save t-axis floats into new image
echo 3drefit -saveatr -atrfloat TAXIS_FLOATS "$taxis_floats" $epi.$corrstr+orig
     3drefit -saveatr -atrfloat TAXIS_FLOATS "$taxis_floats" $epi.$corrstr+orig
#copy over and save t-axis offsets into new image registered dataset
echo 3drefit -saveatr -atrfloat TAXIS_OFFSETS "$taxis_offset" $epi.$corrstr+orig
     3drefit -saveatr -atrfloat TAXIS_OFFSETS "$taxis_offset" $epi.$corrstr+orig
echo 3dNotes -h "slicewise_moco_inplane.sh $epi+orig" $epi.$corrstr+orig
     3dNotes -h "slicewise_moco_inplane.sh $epi+orig" $epi.$corrstr+orig

#echo 3dcopy $epi.slicemocoxy_afni+orig ../$epi.slicemocoxy_afni+orig -overwrite
#echo 3dcopy $epi.slicemocoxy_afni+orig ../$epi.slicemocoxy_afni+orig -overwrite >> slomoco_history.txt
#     3dcopy $epi.slicemocoxy_afni+orig ../$epi.slicemocoxy_afni+orig -overwrite

echo "End of SLOMOCO script" >> slomoco_history.txt
echo "`date`" >> slomoco_history.txt
echo "" >> slomoco_history.txt


