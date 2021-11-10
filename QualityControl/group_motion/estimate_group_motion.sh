#!/bin/bash
# Clear terminal screen
#printf "\033c"

SUB=$1
CONFIG=$2 #/home/fsluser/Documents/rs_proc/conf_files/rep_impact_belgium.conf

# PARSE Configuration file
ABIN=$(awk -F\=  '/^ABIN/{print $2}' "${CONFIG}" | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')  
MNI_REF=$(awk -F \= '/^\<MNI_REF\>/{print $2}' "${CONFIG}" | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')  
MNI_REF_2mm=$(awk -F\=  '/^\<MNI_REF_2mm\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')  
TMPDIR=$(awk -F\=  '/^\<TMPDIR\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')  
OUTDIR=$(awk -F\=  '/^\<OUTDIR\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')  
FINALDIR=$(awk -F\=  '/^\<FINALDIR\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')  
DO_COPY=$(awk -F\=  '/^\<DO_COPY\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')  
DO_REG=$(awk -F\=  '/^\<DO_REG\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')  
DO_REG2=$(awk -F\=  '/^\<DO_REG2\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')  
DO_ICA=$(awk -F\=  '/^\<DO_ICA\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')  
DO_NORM=$(awk -F\=  '/^\<DO_NORM\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')  
DO_QA=$(awk -F\=  '/^\<DO_QA\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')  
DO_QA_NATIVE=$(awk -F\=  '/^\<DO_QA_NATIVE\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')  
DO_CLEAN=$(awk -F\=  '/^\<DO_CLEAN\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')  
TR=$(awk -F\=  '/^\<TR\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')  
ACQ_TYPE=$(awk -F\=  '/^\<ACQ_TYPE\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
USE_BBR=$(awk -F\=  '/^\<USE_BBR\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
DO_SLC=$(awk -F\=  '/^\<DO_SLC\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
DO_FMAP=$(awk -F\=  '/^\<DO_FMAP\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
DO_PEST=$(awk -F\=  '/^\<DO_PEST\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
DO_ANAT_FMAP=$(awk -F\=  '/^\<DO_ANAT_FMAP\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
CORRECT_ANAT_FMAP_NATIVE=$(awk -F\=  '/^\<CORRECT_ANAT_FMAP_NATIVE\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
DO_DEOBL=$(awk -F\=  '/^\<DO_DEOBL\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
DO_DPK=$(awk -F\=  '/^\<DO_DPK\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
DO_ING=$(awk -F\=  '/^\<DO_ING\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
GLOB_VAL=$(awk -F\=  '/^\<GLOB_VAL\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
BREXT_TYPE=$(awk -F\=  '/^\<BREXT_TYPE\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
DO_SMOOTH=$(awk -F\=  '/^\<DO_SMOOTH\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
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




print_debug()
{
	MSG=$1
	MSG_CMD=$2
	HOUR=`date +%T`
	DAY=`date +%F`
	# This makes easier to find the message in log files
	echo '//////////////////////////////////////////////////////////////////////////////////'
	echo ''

	echo "[DEBUG] (${DAY} - ${HOUR} ) ${MSG}"
	echo "[DEBUG] $MSG_CMD"

	echo ''
	echo '//////////////////////////////////////////////////////////////////////////////////'
}


# =================================================
#
#		PREPROCESSING START
#
# =================================================


# Call the project specific script
source ${IMPORT_SCRIPT} ${SUB} ${OUTDIR} 

if [ ! -d "$OUTDIR" ]; then
	echo "${SUB} does not exist. Skipping."
	exit 1
fi


# ------------------------------------------------------------------------

VOL2MED_COR=`3dTqual -automask ${OUTDIR}/func_data.nii`
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


3dvolreg -linear -prefix ${OUTDIR}/r_func_data.nii -maxite 5 -base ${REF_VOL} -1Dfile ${OUTDIR}/motion_estimate.par  ${OUTDIR}/func_data.nii


python3.6 ~/Documents/rs_proc/QC_funcs/calc_motion.py -in ${OUTDIR} -out ${OUTDIR} -m motion_estimate.par -f r_func_data.nii


rm ${OUTDIR}/r_func_data.nii
rm ${OUTDIR}/func_data.nii


cp ${OUTDIR}/motion_summary.txt ${FINALDIR}/${SUB}_mot.txt
rm -r ${OUTDIR}

