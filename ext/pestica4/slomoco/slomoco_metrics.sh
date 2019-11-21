#! /bin/bash

function Usage () {
  cat <<EOF
   slomoco_metrics (distributed with PESTICA v2.1)

 Usage:  slomoco_metrics.sh -d <epi_filename> -s <slice_timing>
 	     -d = dataset: <epi_filename> is the file prefix of
	     the 3D+time dataset to use and correct.
               Note: this script will detect suffix for epi_filename
	     -s = slice_timing integer (default if not specified is = 1)
	          1=alternating interleaved ascending Siemens order (odd/even for odd#, even/odd for even#)
		  2=sequential ascending
	          3=alternating interleaved ascending GE/Phillips order (odd/even)
	     -f = filter width (default if not specified is = 1)
	          3 or 5 width Savitsky-Golay filter widths are supported

EOF
  exit 1
}

# default slice acq order is 1
SLOMOCO_SLICE_ORDER=1
FILTER_WIDTH=1
while getopts hd:s:f: opt; do
  case $opt in
    h)
       Usage
       exit 1
       ;;
    d) # base 3D+time EPI dataset to use to perform corrections
       epi=$OPTARG
       ;;
    s) # slice timing integer, interleaved asc or sequential ascending supported
       SLOMOCO_SLICE_ORDER=$OPTARG
       ;;
    f) # Savitsky-Golay filter width. Only 1,3, or 5 are supported, default is 1 (no filter)
       FILTER_WIDTH=$OPTARG
       ;;
    :)
      echo "option requires input"
      exit 1
      ;;
  esac
done

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
epi_pestica="$epi"_pestica
if [ ! -d $epi_pestica ] ; then
 echo ""
 echo "* PESTICA directory does not exist, SLOMOCO has not been run yet, exiting"
 echo ""
 exit 1
fi

echo "*****   Using copied input in $epi_pestica/$epi+orig.HEAD as input timeseries"
# always copy input file into PESTICA subdirectory in AFNI format
if [ ! -f $epi_pestica/$epi+orig.HEAD ] ; then
  echo "Copying: 3dcopy $epi$suffix $epi_pestica/$epi"
  3dcopy $epi$suffix $epi_pestica/$epi
  echo "3dcopy $epi$suffix $epi_pestica/$epi" >> $epi_pestica/pestica_history.txt
  echo "" >> $epi_pestica/pestica_history.txt
fi

# write command line and SLCMOCO_DIR to history file
echo "`date`" >> $epi_pestica/pestica_history.txt
echo "PESTICA_afni-v2 command line: `basename $fullcommand` $*" >> $epi_pestica/pestica_history.txt
echo "PESTICA env:
`env | grep PESTICA`" >> $epi_pestica/pestica_history.txt
echo "" >> $epi_pestica/pestica_history.txt

# do ALL work inside the pestica/ subdirectory
homedir=`pwd`
cd $epi_pestica

matlab $MATLABLINE <<<"addpath $SLCMOCO_DIR; qa_slomoco('$epi.slicemocoxy_afni+orig','tempslmoco_volslc_alg_vol_$epi.slicemocoxy_afni/motion.wholevol_zt','tempslmocoxy_afni_$epi',$SLOMOCO_SLICE_ORDER,$FILTER_WIDTH); exit;"

cd $homedir

echo "End of PESTICA script" >> $epi_pestica/pestica_history.txt
echo "`date`" >> $epi_pestica/pestica_history.txt
echo "" >> $epi_pestica/pestica_history.txt


