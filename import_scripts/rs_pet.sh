#!/bin/bash

# Copies necessary files to processing [OUTDIR] directory
# !!!IMPORTANT!!!
#
# This is what OUTDIR is expected to contain after importing
#
# t1.nii  	 -> Highres anatomical file
# func_data.nii  -> 4D EPI file
#
# [OPTIONAL]
#
# rev_phase.nii  -> 4D EPI acquired with reverse phase [for fieldmap correction]
SUB=$1
OUTDIR=$2
DATADIR=/media/thiago/EXTRALINUX/fmri_proc/DATA/rsfMRI/$SUB
SPMDIR=/media/thiago/EXTRALINUX/fmri_proc/matlab/toolbox/spm12
SPMDIR2=/media/thiago/EXTRALINUX/fmri_proc/matlab/toolbox/spm12/compat
SPMDIR3=/media/thiago/EXTRALINUX/fmri_proc/matlab/toolbox/spm12/toolbox/OldNorm

if [ ! -d "$DATADIR" ]; then
	echo "Cannot find ${SUB}. Data import failed."
else
	if [ ! -d "${OUTDIR}" ]; then
		mkdir -p -m 777 ${OUTDIR}
	else
	{
		rm ${OUTDIR}/*
	} &> /dev/null
  fi


  IDX=0
	#for F in "${DATADIR}/orig_files/*resting*.nii"
  for F in $(find ${DATADIR}/orig_files/ -type f);
	do

    if [ "$IDX" -lt "10" ]; then
  		C=00${IDX}
  	elif [ "$IDX" -lt "100" ]; then
      C=0${IDX}
    else
  		C=${IDX}
  	fi

    C=$(echo $F| rev | cut -d'_' -f 1 | rev )
    SUFF=$(echo $C| cut -d'.' -f 2  )
    C=$(echo $C| cut -d'.' -f 1  )


    if [ $SUFF == 'nii' ]; then
		  cp $F ${OUTDIR}/tmp_func_data_${C}.nii

      IDX=$((IDX+1))
    fi

	done

  matlab "-nodesktop -nosplash " <<<"addpath('$SPMDIR'); addpath('$SPMDIR2'); addpath('$SPMDIR3'); cd './matlab'; p_reorient('${OUTDIR}/tmp_func_data_', '%03d.nii', 320); exit;" #320

  3dTcat -prefix ${OUTDIR}/func_data.nii -TR 1.7 ${OUTDIR}/tmp_*.nii

  rm -f ${OUTDIR}/tmp_*.nii

  # SET Correct orientation
  #fslorient -copyqform2sform  ${OUTDIR}/func_data.nii
  #3dresample -prefix ${OUTDIR}/tmp.nii -orient ras  -input ${OUTDIR}/func_data.nii
  #mv ${OUTDIR}/tmp.nii ${OUTDIR}/func_data.nii

  # COPY T1 File
	for F in "${DATADIR}/T1/*MR*.nii"
	do
		cp $F ${OUTDIR}/t1.nii

    # SET Correct orientation
    matlab "-nodesktop -nosplash " <<<"addpath('$SPMDIR'); addpath('$SPMDIR2'); addpath('$SPMDIR3'); cd './matlab'; p_reorient('${OUTDIR}/t1.nii'); exit;"
    #mv ${OUTDIR}/temp.nii ${OUTDIR}/t1.nii
    #fslorient -copyqform2sform  ${OUTDIR}/t1.nii
    #3dresample -prefix ${OUTDIR}/tmp.nii -orient ras  -input ${OUTDIR}/t1.nii
    #mv ${OUTDIR}/tmp.nii ${OUTDIR}/t1.nii
	done

  INFO=`fslinfo ${OUTDIR}/t1.nii`
  SX=`sed -n 7p <<< "$INFO"`
  SX=`echo $SX | awk -F ' ' '{print $2}'`

  SY=`sed -n 8p <<< "$INFO"`
  SY=`echo $SY | awk -F ' ' '{print $2}'`

  SZ=`sed -n 9p <<< "$INFO"`
  SZ=`echo $SZ | awk -F ' ' '{print $2}'`

  matlab "-nodesktop -nosplash " <<<"addpath('$SPMDIR'); addpath('$SPMDIR2'); addpath('$SPMDIR3'); cd './matlab'; reslice_images('${OUTDIR}/t1.nii',{'${OUTDIR}/t1.nii'}, 1, 0, [${SX} ${SY} ${SZ}]); exit;"
  3dAutobox -prefix ${OUTDIR}/t1_rb.nii -npad 5 -input ${OUTDIR}/t1_r.nii
  mv ${OUTDIR}/t1_rb.nii ${OUTDIR}/t1.nii
  rm -f ${OUTDIR}/t1_*.nii



  INFO=`fslinfo ${OUTDIR}/func_data.nii`
  SX=`sed -n 7p <<< "$INFO"`
  SX=`echo $SX | awk -F ' ' '{print $2}'`

  SY=`sed -n 8p <<< "$INFO"`
  SY=`echo $SY | awk -F ' ' '{print $2}'`

  SZ=`sed -n 9p <<< "$INFO"`
  SZ=`echo $SZ | awk -F ' ' '{print $2}'`

  matlab "-nodesktop -nosplash " <<<"addpath('$SPMDIR'); addpath('$SPMDIR2'); addpath('$SPMDIR3');  cd  './matlab'; reslice_images('${OUTDIR}/func_data.nii',{'${OUTDIR}/func_data.nii'}, 1, 0, [${SX} ${SY} ${SZ}]); exit;"
  3dAutobox -prefix ${OUTDIR}/func_data_rb.nii -npad 5 -input ${OUTDIR}/func_data_r.nii
  mv ${OUTDIR}/func_data_rb.nii ${OUTDIR}/func_data.nii
  rm -f ${OUTDIR}/func_data_*.nii
fi
