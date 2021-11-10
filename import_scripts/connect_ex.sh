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

print_timestamp_prefix()
{
  F=$1
  IFS='/'
  FNAME=''
  read -ra ADDR <<< "$F"
  for TOK in "${ADDR[@]}"; do
      FNAME=${TOK}
  done

  # 2. Get the timestamp
  IFS='_'
  read -ra ADDR <<< "$FNAME"
  for TOK in "${ADDR[@]}"; do
      RS_TIMESTAMP=${TOK}
      break
  done
  IFS=' ' # reset to default value after usage
  echo "${RS_TIMESTAMP}"
}

SUB=$1
OUTDIR=$2
DATADIR="/media/u0101486/Seagate Backup Plus Drive/ConnectExData/Thiago/$SUB/"

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



  for F in "${DATADIR}"/fMRI/*.nii.gz
  do
    # 1. Get the filename
    RS_TIMESTAMP=`print_timestamp_prefix "${F}"`

  	cp "$F" ${OUTDIR}/func_data.nii.gz
  	gunzip -f ${OUTDIR}/func_data.nii.gz
  done


	for F in "${DATADIR}"/fMRI_rev/*rsf*.nii.gz
	do
		cp "$F" ${OUTDIR}/rev_phase.nii.gz
		gunzip -f ${OUTDIR}/rev_phase.nii.gz
	done


	for F in "${DATADIR}"/T1/*MPRAGE*.nii.gz
	do
		cp "$F" ${OUTDIR}/t1.nii.gz
		gunzip -f ${OUTDIR}/t1.nii.gz
	done



  CURRENT_LARGE=0
  LARGEST_LOG=0

  # Find out which log is the largest
  for F in "${DATADIR}"/${SUB}_phys/*.log
  do
    LOG_TIMESTAMP=`stat -t -c %s "${F}"`

    if [ "${LOG_TIMESTAMP}" -gt "${CURRENT_LARGE}" ]; then
      CURRENT_LARGE=${LOG_TIMESTAMP}
      LARGEST_LOG=${F}
    fi
  done


  # Copy the largest log file
  # Generally, this would correspond to the functional data
  # (this certainly is true for the present dataset)
  # The following line will raise an error if there are no log files,
  # though this would be a different problem altogether
  cp -f "${LARGEST_LOG}" ${OUTDIR}/physio.log
fi
