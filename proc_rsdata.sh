#!/bin/bash
# Clear terminal screen
#printf "\033c"

SUB=$1
CONFIG=$2 #/home/fsluser/Documents/rs_proc/conf_files/rep_impact_belgium.conf
BLOCK=$3

QCDIR='/home/fsluser/Documents/rs_proc/QualityControl/'

# PARSE Configuration file
ABIN=$(awk -F\=  '/^ABIN/{print $2}' "${CONFIG}" | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')  
FIXBIN=$(awk -F\=  '/^FIXBIN/{print $2}' "${CONFIG}" | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')  
MNI_REF=$(awk -F \= '/^\<MNI_REF\>/{print $2}' "${CONFIG}" | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')  
MNI_REF_2mm=$(awk -F\=  '/^\<MNI_REF_2mm\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')  
TMPDIR=$(awk -F\=  '/^\<TMPDIR\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')  
OUTDIR=$(awk -F\=  '/^\<OUTDIR\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')  
FINALDIR=$(awk -F\=  '/^\<FINALDIR\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')  
DO_COPY=$(awk -F\=  '/^\<DO_COPY\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')  
MOCO=$(awk -F\=  '/^\<MOCO\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
DO_REG=$(awk -F\=  '/^\<DO_REG\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')  
DO_REG2=$(awk -F\=  '/^\<DO_REG2\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')  
DO_ICA=$(awk -F\=  '/^\<DO_ICA\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')  
DO_NORM=$(awk -F\=  '/^\<DO_NORM\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')  
DO_QA=$(awk -F\=  '/^\<DO_QA\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')  
DO_QA_NATIVE=$(awk -F\=  '/^\<DO_QA_NATIVE\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')  
DO_CLEAN=$(awk -F\=  '/^\<DO_CLEAN\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')  
DO_FUNC_BIAS=$(awk -F\=  '/^\<DO_FUNC_BIAS\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')  
DO_ATLAS_NAT=$(awk -F\=  '/^\<DO_ATLAS_NAT\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')  
TR=$(awk -F\=  '/^\<TR\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')  
ACQ_TYPE=$(awk -F\=  '/^\<ACQ_TYPE\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
USE_BBR=$(awk -F\=  '/^\<USE_BBR\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
DO_SLC=$(awk -F\=  '/^\<DO_SLC\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
DO_FMAP=$(awk -F\=  '/^\<DO_FMAP\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
EES=$(awk -F\=  '/^\<EES\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
DO_PEST=$(awk -F\=  '/^\<DO_PEST\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
DO_ANAT_FMAP=$(awk -F\=  '/^\<DO_ANAT_FMAP\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
CORRECT_ANAT_FMAP_NATIVE=$(awk -F\=  '/^\<CORRECT_ANAT_FMAP_NATIVE\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
DO_DEOBL=$(awk -F\=  '/^\<DO_DEOBL\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
DO_DPK=$(awk -F\=  '/^\<DO_DPK\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
DO_ING=$(awk -F\=  '/^\<DO_ING\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
GLOB_VAL=$(awk -F\=  '/^\<GLOB_VAL\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')

FIX_CLASSIFIER=$(awk -F\=  '/^\<FIX_CLASSIFIER\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
FIX_CL_LABEL=$(awk -F\=  '/^\<FIX_CL_LABEL\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
FIX_THR=$(awk -F\=  '/^\<FIX_THR\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')


BREXT_TYPE=$(awk -F\=  '/^\<BREXT_TYPE\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
#DEPRECATED DO_SMOOTH=$(awk -F\=  '/^\<DO_SMOOTH\>/{print $2}' "${CONFIG}"  | cut -d'#' -f1 |  sed -e 's/[[:space:]]*$//')
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




# In case data has multiple blocks (multiple resting state blocks or perhaps task)
if test ! -z "${BLOCK}"; then
	FINALDIR=${FINALDIR}_${BLOCK}
	OUTDIR=${OUTDIR}_${BLOCK}
fi

#rm -f -r ${OUTDIR}/dicer


#gunzip -f ${OUTDIR}/t1_frn_seg.nii.gz

#mv ${OUTDIR}/t1_frn_seg.nii ${OUTDIR}/t1_segt.nii

#rm -f ${OUTDIR}/proc_data_native_3.nii
#3dresample -dxyz 3 3 3 -prefix ${OUTDIR}/proc_data_native_3.nii -input ${OUTDIR}/proc_data_native.nii
#
#${ABIN}/antsApplyTransforms \
#	-i ${OUTDIR}/t1_frn.nii \
#	--float \
#	-r ${OUTDIR}/proc_data_native_3.nii \
#	-o ${OUTDIR}/t1_nat.nii \
#	-t [${OUTDIR}/func2anat.mat,1] \
#	-n Linear
#
#CURDIR=`pwd`
#cd /home/fsluser/Documents/DiCER-master/
#mkdir ${OUTDIR}/dicer
#
#
#source DiCER_lightweight.sh -i${OUTDIR}/proc_data_native_3.nii -a ${OUTDIR}/t1_nat.nii -w ${OUTDIR}/dicer -s SUBJECT_1  -p 2 -d 
#cd ${CURDIR}
#exit 1

print_debug()
{
	MSG=$1
	MSG_CMD=$2
	HOUR=`date +%T`
	DAY=`date +%F`
	# This makes easier to find the message in log files
	echo ''
	echo ''
	echo '//////////////////////////////////////////////////////////////////////////////////'
	echo ''

	echo "[DEBUG] (${DAY} - ${HOUR} ) ${MSG}"
	echo "$MSG_CMD"

	echo ''
	echo '//////////////////////////////////////////////////////////////////////////////////'
	echo ''
	echo ''
}


log_command_div()
{
	HOUR=`date +%T`
	DAY=`date +%F`
	
	echo ''
	echo ''
	echo '*****************************************************************************************'
	echo "(${DAY} - ${HOUR} ) COMMAND DONE"
	echo '*****************************************************************************************'
	echo ''
	echo ''

}




# =================================================
#
#		PREPROCESSING START
#
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



			if (( $TN < 10 )); then
				FI=00$TN
			elif (( $TN < 100 )); then
				FI=0$TN
			else
				FI=$TN	
			fi

			{
			3dAutomask -apply_prefix ${OUTDIR}/tmp/func_data_m.${FI}.nii -prefix ${OUTDIR}/tmp/mask.${FI}.nii ${OUTDIR}/tmp/func_data.${FI}.nii
			${ABIN}/DenoiseImage -s 1 -x ${OUTDIR}/tmp/mask.${FI}.nii -n Rician -v 1 -i ${OUTDIR}/tmp/func_data_m.${FI}.nii  -o [ ${OUTDIR}/tmp/func_data_n.${FI}.nii ]

			} &> /dev/null

		}
		export -f skull_strip


		print_debug 'Performing volume-wise automasking'

		NVOLS=`fslnvols ${OUTDIR}/func_data.nii`
		NZ=`3dinfo -nk ${OUTDIR}/func_data.nii`

		print_debug "NVOLS = ${NVOLS}"
		python3.6 ~/Documents/rs_proc/util/write_slice_timing.py -tr ${TR} -nsl ${NZ} -acq ${ACQ_TYPE} -out ${OUTDIR}/slice_acq.txt

		# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		#
		#  Bias field correction
		#
		#  It is highly unlikely that this step is necessary, but in some scanners
		#  the smooth difference in intensities might affect motion correction
		#  I still need to test this, though
		#
		# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		if [ "${DO_FUNC_BIAS}" -eq "1" ]; then
			log_command_div
			PREF='b'
			mkdir ${OUTDIR}/tmp
			3dTsplit4D -prefix ${OUTDIR}/tmp/func_data.nii -keep_datum ${OUTDIR}/func_data.nii
			parallel -j4 --line-buffer estimate_func_bias ::: $(seq 0 ${NVOLS}) ::: ${OUTDIR} ::: ${ABIN} 

			3dTcat -prefix ${OUTDIR}/func_biasfield.nii -tr ${TR} ${OUTDIR}/tmp/biasfield_*.nii

			fslmaths ${OUTDIR}/func_biasfield.nii -Tmean ${OUTDIR}/mean_biasfield.nii
			fslmaths ${OUTDIR}/func_data.nii -div ${OUTDIR}/mean_biasfield.nii ${OUTDIR}/${PREF}func_data.nii
			gunzip -f ${OUTDIR}/${PREF}func_data.nii.gz

			rm -r ${OUTDIR}/tmp
			python3.6 ${QCDIR}/save_image_diff.py -o ${OUTDIR} -i ${OUTDIR} -a func_data.nii -b ${PREF}func_data.nii -type mean -msg1 'With Bias' -msg2 'Without Bias'
			log_command_div
		else 
			PREF=''
		fi
		
		# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		# Skull stripping
		# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		log_command_div
		mkdir ${OUTDIR}/tmp
		3dTsplit4D -prefix ${OUTDIR}/tmp/func_data.nii -keep_datum ${OUTDIR}/${PREF}func_data.nii
		parallel -j4 --line-buffer skull_strip ::: $(seq 0 ${NVOLS}) ::: ${OUTDIR} ::: ${ABIN} 

		PREF=m${PREF}
		3dTcat -prefix ${OUTDIR}/${PREF}func_data.nii -tr ${TR} ${OUTDIR}/tmp/func_data_n*.nii

		rm -r ${OUTDIR}/tmp
		log_command_div


		# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#		if [ "$DO_DPK" -eq "1" ]; then
#			print_debug 'Despiking [-localedit -NEW]'
#			3dDespike -nomask -NEW -localedit -prefix ${OUTDIR}/d${PREF}func_data.nii ${OUTDIR}/${PREF}func_data.nii 
#			python3.6 ${QCDIR}/save_image_diff.py -o ${OUTDIR} -i ${OUTDIR} -a ${PREF}func_data.nii -b d${PREF}func_data.nii -msg1 'Before Despike' -msg2 'After Despike'
#			PREF=d${PREF}
#			log_command_div
#		fi



		



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
#			3dvolreg -heptic -prefix ${OUTDIR}/r${PREF}func_data.nii -base ${REF_VOL} -rot_thresh 0.01 -delta 3 \
#				-x_thresh 0.01 -zpad 5 -maxite 75 -1Dfile ${OUTDIR}/motion_estimate.par -maxdisp1D ${OUTDIR}/maximum_disp.1d \
#				 ${OUTDIR}/${PREF}func_data.nii

			3dvolreg -heptic -prefix ${OUTDIR}/tmp_1.nii -base ${REF_VOL} -rot_thresh 0.03 -delta 5 \
				-x_thresh 0.03 -zpad 5 -maxite 15 -1Dfile ${OUTDIR}/motion_estimate.par -maxdisp1D ${OUTDIR}/maximum_disp.1d \
				 ${OUTDIR}/${PREF}func_data.nii

			3dAutomask -prefix ${OUTDIR}/mask.nii ${OUTDIR}/tmp_1.nii
			matlab  "-nodesktop -nosplash -softwareopengl " <<<"despike_data('${OUTDIR}', 'tmp_1.nii', 'tmp_2.nii', 'motion_estimate.par', 'mask.nii' ); exit;"
			3dDespike -nomask -NEW -localedit -prefix ${OUTDIR}/tmp_3.nii ${OUTDIR}/tmp_2.nii 



			3dvolreg -heptic -prefix ${OUTDIR}/r${PREF}func_data.nii -base 0 -rot_thresh 0.01 -delta 0.5 \
				-x_thresh 0.01 -zpad 3 -maxite 60 -1Dfile ${OUTDIR}/motion_estimate_2.par -maxdisp1D ${OUTDIR}/maximum_disp_2.1d \
				 ${OUTDIR}/tmp_3.nii

			rm -f ${OUTDIR}/tmp_1.nii
			rm -f ${OUTDIR}/tmp_2.nii
			rm -f ${OUTDIR}/tmp_3.nii


			PREF=r${PREF}



			1dplot -thick -volreg -png ${OUTDIR}/motion_estimate.png -one ${OUTDIR}/motion_estimate.par
			1dplot -thick -png ${OUTDIR}/maximum_disp.png -one ${OUTDIR}/maximum_disp.1d
			1dplot -thick -png ${OUTDIR}/maximum_disp_delt.png -one ${OUTDIR}/maximum_disp.1d_delt

			1dplot -thick -volreg -png ${OUTDIR}/motion_estimate_2.png -one ${OUTDIR}/motion_estimate_2.par
			1dplot -thick -png ${OUTDIR}/maximum_disp_2.png -one ${OUTDIR}/maximum_disp_2.1d
			1dplot -thick -png ${OUTDIR}/maximum_disp_delt_2.png -one ${OUTDIR}/maximum_disp_2.1d_delt

			log_command_div

		fi

		if [ "$MOCO" == "slomoco" ]; then
			print_debug 'DOING SLOMOCO'
			source ~/Documents/rs_proc/align_pestica.sh ${OUTDIR} ${PREF}func_data ${TR} 1



			#TODO Remove this line [only here for comparison]
			3dvolreg -heptic -prefix ${OUTDIR}/r2${PREF}func_data.nii -base 0 -rot_thresh 0.01 -delta 2 \
				-x_thresh 0.01 -zpad 10 -maxite 60 -1Dfile ${OUTDIR}/motion_estimate.par -maxdisp1D ${OUTDIR}/maximum_disp.1d \
				 ${OUTDIR}/p${PREF}func_data.nii

		
			3drefit -deoblique ${OUTDIR}/r2${PREF}func_data.nii
			1dplot -thick -volreg -png ${OUTDIR}/motion_estimate.png -one ${OUTDIR}/motion_estimate.par
			1dplot -thick -png ${OUTDIR}/maximum_disp.png -one ${OUTDIR}/maximum_disp.1d
			1dplot -thick -png ${OUTDIR}/maximum_disp_delt.png -one ${OUTDIR}/maximum_disp.1d_delt


			python3.6 ${QCDIR}/save_image_diff.py -o ${OUTDIR} -i ${OUTDIR} -a r2${PREF}func_data.nii -b rp${PREF}func_data.nii -msg1 'PESTICA+3dvorleg' -msg2 'PESTICA+SLOMOCO'
			cp ${OUTDIR}/slomoco4/*.txt ${OUTDIR}/
			cp ${OUTDIR}/slomoco4/*.1D ${OUTDIR}/
		
			PREF=rp${PREF}
			log_command_div
		fi
		


		# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		if [ "$DO_PEST" -eq "1" ] && [ ! "$MOCO" == "slomoco" ]; then
			print_debug 'PESTICA4 cardiac/respiratory effects correction'
			# 0 here indicates that SLOMOCO procedure should not be done
			source ~/Documents/rs_proc/align_pestica.sh ${OUTDIR} ${PREF}func_data ${TR} 0
			python3.6 ${QCDIR}/save_image_diff.py -o ${OUTDIR} -i ${OUTDIR} -a ${PREF}func_data.nii -b p${PREF}func_data.nii -msg1 'Before PESTICA' -msg2 'After PESTICA'
			PREF=p${PREF}
			log_command_div
		fi

		# Create native space brain mask
		3dAutomask -prefix ${OUTDIR}/nat_mask.nii ${OUTDIR}/${PREF}func_data.nii


		# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		if [ "$DO_SLC" -eq "1" ]; then
			print_debug 'Slice timing correction'
			# Slice timing correction  [Parker et al. 2017 -- 10.1016/j.media.2016.08.006]
			CF=`echo "(1/${TR})/2" | bc -l`

			if (( $(echo "$TR >= 2.5" | bc -l) )); then
				# Used to preserve the full power spectrum
				CF=0.21
			fi

			# Correct slice acquisition to the middle of the volume acquisition
			# To correct to the first slice, set REF_TIME to 0
			# TODO [10.09.19] Pass this to the parameter section of the file
			REF_TIME=`echo "${TR}/2" | bc -l`


			print_debug 'Slice timing correction' "filtershift --in=${OUTDIR}/${PREF}func_data.nii --itl=${ACQ_TYPE}  --cf=${CF} --rt=${REF_TIME} --TR=${TR} --out=${OUTDIR}/a${PREF}func_data.nii"

			filtershift --in=${OUTDIR}/${PREF}func_data.nii --itl=${ACQ_TYPE}  --cf=${CF} --rt=${REF_TIME} --TR=${TR} --out=${OUTDIR}/a${PREF}func_data.nii
			gunzip -f ${OUTDIR}/a${PREF}func_data.nii.gz

			PREF=a${PREF}
			log_command_div
		fi


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

			# Fieldmap estimation [ONLY estimation]
			print_debug 'TOPUP command' "topup --imain=${OUTDIR}/topup_data.nii --datain=${TMPDIR}/topupfield.txt --out=${OUTDIR}/topup_results --fout=${OUTDIR}/fieldmap --iout=${OUTDIR}/w_topup_data --estmov=0,0,0  --regmod=membrane_energy  --minmet=1,1,1 --verbose --warpres=8,6,4 --miter=32,12,4 --subsamp=4,2,1 --fwhm=6,4,2"

			topup --imain=${OUTDIR}/topup_data.nii --datain=${TMPDIR}/topupfield.txt --out=${OUTDIR}/topup_results --fout=${OUTDIR}/fieldmap --iout=${OUTDIR}/w_topup_data \
				--estmov=0,0,0  --regmod=membrane_energy  --minmet=1,1,1 --verbose --warpres=8,6,4 --miter=32,12,4 --subsamp=2,2,1 --fwhm=6,4,2
			
			applytopup --imain=${OUTDIR}/${PREF}func_data.nii --datain=${TMPDIR}/topupfield.txt --topup=${OUTDIR}/topup_results --inindex=1 --method=jac --interp=spline --out=${OUTDIR}/u${PREF}func_data
			gunzip -f ${OUTDIR}/u${PREF}func_data.nii.gz
			3dAutomask -apply_prefix ${OUTDIR}/tmp.nii ${OUTDIR}/u${PREF}func_data.nii
			rm ${OUTDIR}/u${PREF}func_data.nii
			mv ${OUTDIR}/tmp.nii ${OUTDIR}/u${PREF}func_data.nii
			python3.6 ${QCDIR}/save_image_diff.py -o ${OUTDIR} -i ${OUTDIR} -a ${PREF}func_data.nii -b u${PREF}func_data.nii -type mean -msg1 'Before FMAP' -msg2 'After FMAP'
			PREF=u${PREF}
			log_command_div

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
			python3.6 ${QCDIR}/save_image_diff.py -o ${OUTDIR} -i ${OUTDIR} -a ${PREF}func_data.nii -b u${PREF}func_data.nii -type mean -msg1 'Before FMAP' -msg2 'After FMAP'

			PREF=u${PREF}
			log_command_div
		fi


		# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		if [ "$DO_DEOBL" -eq "1" ]; then
			print_debug 'Deobliquing volumes'
			#3dWarp -deoblique -newgrid ${DEOBL_VOX} -NN -prefix ${OUTDIR}/w${PREF}func_data.nii ${OUTDIR}/${PREF}func_data.nii
			cp ${OUTDIR}/${PREF}func_data.nii ${OUTDIR}/w${PREF}func_data.nii
			3drefit -deoblique ${OUTDIR}/w${PREF}func_data.nii
			PREF=w${PREF}
			log_command_div
		fi



		


		# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		

		if [ "$DO_ING" -eq "1" ]; then
			print_debug 'Global intensity normalization' 'fslmaths ${OUTDIR}/${PREF}func_data -ing ${GLOB_VAL} ${OUTDIR}/i${PREF}func_data'
			fslmaths ${OUTDIR}/${PREF}func_data -ing ${GLOB_VAL} ${OUTDIR}/i${PREF}func_data
			gunzip -f ${OUTDIR}/i${PREF}func_data.nii.gz
			PREF=i${PREF}
		fi


		# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		# Update brain mask
		rm ${OUTDIR}/nat_mask.nii
		3dAutomask -prefix ${OUTDIR}/nat_mask.nii -apply_prefix ${OUTDIR}/proc_data_native.nii ${OUTDIR}/${PREF}func_data.nii

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
if [ "$DO_REG2" -eq "1" ]; then
	printf "[$SUB] Functional <-> Anat <-> MNI Registration ... "
	START=$(date -u +%s.%N)


	{

	# Bias field correction and Gaussian noise reduction
	# TODO Pass to parameter section
	# Necessary files must be copied in case this section does not run (e.g. multiple blocks of functional data for the same participant)
	DO_ANAT=1

	if [ "$DO_ANAT" -eq "1" ]; then
		#3dWarp -deoblique -newgrid ${DEOBL_VOX_ANAT} -NN -prefix ${OUTDIR}/t1_f.nii ${OUTDIR}/t1.nii
		cp ${OUTDIR}/t1.nii ${OUTDIR}/t1_f.nii
		if [ "$DO_DEOBL" -eq "1" ]; then
			3drefit -deoblique ${OUTDIR}/t1_f.nii
		fi

		${ABIN}/N4BiasFieldCorrection -i ${OUTDIR}/t1_f.nii -s 4  -o [ ${OUTDIR}/t1_fb.nii,${OUTDIR}/biasfield.nii ]

		
		if [ "$BREXT_TYPE" -eq "1" ]; then
			python3.6 ~/Documents/rs_proc/extract_brain.py -i "${OUTDIR}/t1_fb.nii" -o "${OUTDIR}/AnatMask.nii"
		fi


		if [ "$BREXT_TYPE" -eq "2" ]; then
			#TODO
			echo 'TODO'
		fi

		3dcalc -a ${OUTDIR}/t1_fb.nii -b ${OUTDIR}/AnatMask.nii -expr 'a*b' -prefix ${OUTDIR}/t1_brain.nii
		${ABIN}/DenoiseImage -i ${OUTDIR}/t1_brain.nii  -o [ ${OUTDIR}/t1_frn.nii,${OUTDIR}/t1_noise.nii ]

	

		# Anat2MNI. The results of this transform are also used to obtain priors for segmentation
		SRC=${OUTDIR}/t1_frn.nii

		${ABIN}/antsRegistration -d 3 -r [$MNI_REF,$SRC,0] -v 1 \
				-m MI[$MNI_REF,$SRC,1,32] -t translation[0.1] -c [500,5.e-7,20] \
				-s 3vox -f 3 -l 1 -n BSpline \
				-m MI[$MNI_REF,$SRC,1,32,Regular,0.25] -t rigid[0.1] -c [500,5.e-7,20] \
				-s 3vox -f 3 -l 1 -n BSpline \
				-m MI[$MNI_REF,$SRC,1,32,Regular,0.25] -t affine[0.1] -c [500x100x10,5.e-7,10] \
				-s 2x1x0vox -f 3x2x1 -l 1 -n BSpline \
				-m CC[$MNI_REF,$SRC,1,3] -t SyN[0.1,3] -c [25x10,1.e-7,10] \
				-s 1x0vox -f 2x1 -l 1 -n BSpline \
				-o [${OUTDIR}/anat2group,${OUTDIR}/anat2group.nii]


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

	fi

	# Prepare reference functional image
	# Strictly speaking, it is not necessary to denoise, but my testes indicate that this might help in some cases
	# I cannot see a reason it would negatively impact anything
	fslmaths "${OUTDIR}/proc_data_native.nii" -Tmedian -thr 0 "${OUTDIR}/ref_func_b.nii"
	gunzip -f ${OUTDIR}/ref_func_b.nii.gz
	${ABIN}/DenoiseImage -i ${OUTDIR}/ref_func_b.nii  -o [ ${OUTDIR}/ref_func_bf.nii]
	3dAutomask -apply_prefix ${OUTDIR}/ref_func_bfm.nii ${OUTDIR}/ref_func_bf.nii 



	if [ "${USE_BBR}" -eq "1" ]; then

		echo "Using BBR procedure to match func to anat images"
		epi_reg --epi=${OUTDIR}/ref_func_bf.nii --t1brain=${OUTDIR}/t1_frn.nii --t1=${OUTDIR}/t1_f.nii \
			--wmseg=${OUTDIR}/wm_mask.nii -v --out=${OUTDIR}/afunc2anat.nii 
		gunzip -f ${OUTDIR}/afunc2anat.nii.gz

	else
		flirt -in ${OUTDIR}/ref_func_bfm.nii -out ${OUTDIR}/afunc2anat.nii -ref ${OUTDIR}/t1_frn.nii -cost normmi -omat ${OUTDIR}/afunc2anat.mat 
		gunzip -f ${OUTDIR}/afunc2anat.nii
	fi

	~/Documents/itk/c3d/bin/c3d_affine_tool -ref ${OUTDIR}/t1_frn.nii -src ${OUTDIR}/ref_func_bf.nii ${OUTDIR}/afunc2anat.mat -fsl2ras -oitk ${OUTDIR}/func2anat.mat


	
	
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
	

				#TODO Check whether this step is wanted
				${ABIN}/DenoiseImage -s 1 -x ${OUTDIR}/nat_mask.nii -n Rician -v 1 -i ${OUTDIR}/tmp/proc_data_corr_${FI}.nii  -o [ ${OUTDIR}/tmp/proc_data_s.${FI}.nii ]
				mv ${OUTDIR}/tmp/proc_data_s.${FI}.nii ${OUTDIR}/tmp/proc_data.${FI}.nii


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


			
			python3.6 ~/Documents/rs_proc/QC_funcs/save_fmap_comp.py -o ${OUTDIR} -i ${OUTDIR} -a proc_data_native.nii -b proc_data_native_fmap.nii -c t1_func.nii -msg1 'Before FMAP' -msg2 'After FMAP'
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
fi
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------


if [ "$DO_ATLAS_NAT" -eq "1" ]; then

	printf "[$SUB] Registering Atlases to Native Space ... "
	START=$(date -u +%s.%N)
	{
	#TODO [11.10.19] Change REF location to be configurable 
	# Although, this may be unnecessary, as this is related to the atlases and therefore fixed
	REF=/mnt/hgfs/ssd_tmp/atlases/MNI152_T1_1mm_brain.nii
	AAL=/mnt/hgfs/ssd_tmp/atlases/aal/aal2/AAL2.nii
	LOCAL_GLOBAL=/mnt/hgfs/ssd_tmp/atlases/local_global/Schaefer2018_400Parcels_17Networks_order_FSLMNI152_1mm.nii
	SRC=${OUTDIR}/t1_frn.nii

	${ABIN}/antsRegistration -d 3 -r [$MNI_REF,$SRC,0] -v 1 \
			-m MI[$MNI_REF,$SRC,1,32] -t translation[0.1] -c [500,5.e-7,20] \
			-s 3vox -f 3 -l 1 -n BSpline \
			-m MI[$MNI_REF,$SRC,1,32,Regular,0.25] -t rigid[0.1] -c [500,5.e-7,20] \
			-s 3vox -f 3 -l 1 -n BSpline \
			-m MI[$MNI_REF,$SRC,1,32,Regular,0.25] -t affine[0.1] -c [500x100x10,5.e-7,10] \
			-s 2x1x0vox -f 3x2x1 -l 1 -n BSpline \
			-m CC[$MNI_REF,$SRC,1,3] -t SyN[0.1,3] -c [25x10,1.e-7,10] \
			-s 1x0vox -f 2x1 -l 1 -n BSpline \
			-o [${OUTDIR}/anat2atlas,${OUTDIR}/anat2atlas.nii]


	if [ "${CORRECT_ANAT_FMAP_NATIVE}" -eq "0" ] && test -f "${OUTDIR}/func2anat0Warp.nii.gz"; then
		${ABIN}/antsApplyTransforms -v 1 \
			-i ${AAL} \
			--float \
			-r ${OUTDIR}/ref_func_bf.nii \
			-o ${OUTDIR}/aal2_nat_atlas.nii \
			-t [${OUTDIR}/func2anat.mat,1] \
			-t [${OUTDIR}/func2anat0InverseWarp.nii.gz,0] \
			-t [${OUTDIR}/anat2group0GenericAffine.mat,1] \
			-t [${OUTDIR}/anat2atlas1InverseWarp.nii.gz,0] \
			-n NearestNeighbor


		${ABIN}/antsApplyTransforms -v 1 \
			-i ${LOCAL_GLOBAL} \
			--float \
			-r ${OUTDIR}/ref_func_bf.nii \
			-o ${OUTDIR}/local_global400_nat_atlas.nii \
			-t [${OUTDIR}/func2anat.mat,1] \
			-t [${OUTDIR}/func2anat0InverseWarp.nii.gz,0] \
			-t [${OUTDIR}/anat2group0GenericAffine.mat,1] \
			-t [${OUTDIR}/anat2atlas1InverseWarp.nii.gz,0] \
			-n NearestNeighbor
	else
		${ABIN}/antsApplyTransforms -v 1 \
			-i ${LOCAL_GLOBAL} \
			--float \
			-r ${OUTDIR}/ref_func_bf.nii \
			-o ${OUTDIR}/local_global400_nat_atlas.nii \
			-t [${OUTDIR}/func2anat.mat,1] \
			-t [${OUTDIR}/anat2group0GenericAffine.mat,1] \
			-t [${OUTDIR}/anat2atlas1InverseWarp.nii.gz,0] \
			-n NearestNeighbor

		${ABIN}/antsApplyTransforms -v 1 \
			-i ${AAL} \
			--float \
			-r ${OUTDIR}/ref_func_bf.nii \
			-o ${OUTDIR}/aal2_nat_atlas.nii \
			-t [${OUTDIR}/func2anat.mat,1] \
			-t [${OUTDIR}/anat2group0GenericAffine.mat,1] \
			-t [${OUTDIR}/anat2atlas1InverseWarp.nii.gz,0] \
			-n NearestNeighbor
	fi
	} &> ${OUTDIR}/02b_atlas_norm.log

	END=$(date -u +%s.%N)
	DIFF=`echo "( $END - $START )" | bc`
	printf "DONE [%.1f s]\n" $DIFF

fi

# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------


if [ "${REG_MODEL}" != "FIX" ] || [ "$DO_QA" -eq "1" ] || [ "$EXTRACT_NUIS" -eq "1" ]; then

	printf "[$SUB] Extracting Nuisance Regressors ... "
	START=$(date -u +%s.%N)
	{
		
	echo "Using ${REG_MODEL} model for nuisance regression"

	# Warp tissue maps from anatomical to functional
	# If anatomical 'FMAP' correction was not done, only apply affine transform
	if [ "${DO_ANAT_FMAP}" == "1" ]; then

		${ABIN}/antsApplyTransforms -v 1 \
			-i ${OUTDIR}/csf_mask.nii \
			--float \
			-r ${OUTDIR}/ref_func_bf.nii \
			-o ${OUTDIR}/csf_mask_nat.nii \
			-t [${OUTDIR}/func2anat.mat,1] \
			-t [${OUTDIR}/func2anat0InverseWarp.nii.gz,0] \
			-n Linear

		${ABIN}/antsApplyTransforms -v 1 \
			-i ${OUTDIR}/wm_mask.nii \
			--float \
			-r ${OUTDIR}/ref_func_bf.nii \
			-o ${OUTDIR}/wm_mask_nat.nii \
			-t [${OUTDIR}/func2anat.mat,1] \
			-t [${OUTDIR}/func2anat0InverseWarp.nii.gz,0] \
			-n Linear


		${ABIN}/antsApplyTransforms -v 1 \
			-i ${OUTDIR}/gm_mask.nii \
			--float \
			-r ${OUTDIR}/ref_func_bf.nii \
			-o ${OUTDIR}/gm_mask_nat.nii \
			-t [${OUTDIR}/func2anat.mat,1] \
			-t [${OUTDIR}/func2anat0InverseWarp.nii.gz,0] \
			-n Linear
	else

		${ABIN}/antsApplyTransforms -v 1 \
			-i ${OUTDIR}/csf_mask.nii \
			--float \
			-r ${OUTDIR}/ref_func_bf.nii \
			-o ${OUTDIR}/csf_mask_nat.nii \
			-t [${OUTDIR}/func2anat.mat,1] \
			-n Linear

		${ABIN}/antsApplyTransforms -v 1 \
			-i ${OUTDIR}/wm_mask.nii \
			--float \
			-r ${OUTDIR}/ref_func_bf.nii \
			-o ${OUTDIR}/wm_mask_nat.nii \
			-t [${OUTDIR}/func2anat.mat,1] \
			-n Linear

		${ABIN}/antsApplyTransforms -v 1 \
			-i ${OUTDIR}/gm_mask.nii \
			--float \
			-r ${OUTDIR}/ref_func_bf.nii \
			-o ${OUTDIR}/gm_mask_nat.nii \
			-t [${OUTDIR}/func2anat.mat,1] \
			-n Linear
	fi

	# ERODE tissue masks
	fslmaths ${OUTDIR}/csf_mask_nat -thr 0.66 -bin  ${OUTDIR}/csf_mask_nat_ero
	gunzip -f ${OUTDIR}/csf_mask_nat_ero.nii.gz

	fslmaths ${OUTDIR}/wm_mask_nat  -bin -kernel 3D -ero ${OUTDIR}/wm_mask_nat_ero
	gunzip -f ${OUTDIR}/wm_mask_nat_ero.nii.gz

	3dcalc -a ${OUTDIR}/wm_mask_nat_ero.nii -b ${OUTDIR}/csf_mask_nat_ero.nii -expr 'a+b' -prefix ${OUTDIR}/nongm_mask_ero.nii

	fslmeants -i ${OUTDIR}/proc_data_native.nii -m ${OUTDIR}/csf_mask_nat_ero.nii -o ${OUTDIR}/csf_sig.txt
	fslmeants -i ${OUTDIR}/proc_data_native.nii -m ${OUTDIR}/wm_mask_nat_ero.nii -o ${OUTDIR}/wm_sig.txt
	fslmeants -i ${OUTDIR}/proc_data_native.nii -m ${OUTDIR}/nat_mask.nii -o ${OUTDIR}/global_sig.txt


	python3.6 ~/Documents/rs_proc/run_acompcor.py -d ${OUTDIR} -i proc_data_native.nii -n nongm_mask_ero.nii -b nat_mask.nii -t ${TR}

	# Remove headers from compcor results
	sed '1d' ${OUTDIR}/acompcor.txt > ${OUTDIR}/tmp_acompcor.txt
	mv ${OUTDIR}/tmp_acompcor.txt ${OUTDIR}/acompcor.txt

	sed '1d' ${OUTDIR}/tcompcor.txt > ${OUTDIR}/tmp_tcompcor.txt
	mv ${OUTDIR}/tmp_tcompcor.txt ${OUTDIR}/tcompcor.txt

	python3.6 ~/Documents/rs_proc/rp_nuis_calc.py -d ${OUTDIR} -r motion_estimate.par -t ${FD_THR}

	} &> ${OUTDIR}/03_SignalExtraction.log

	END=$(date -u +%s.%N)
	DIFF=`echo "( $END - $START )" | bc`
	printf "DONE [%.1f s]\n" $DIFF

fi 

# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------



# From Smith et al. NeuroImage 2013, HCP

# Preliminary analyses of measures of motion-related artefacts indicate that the ICA-FIX process greatly reduces but does not, in some datasets, totally eliminate motion artefacts 
# that are frame-specific and non-spatially-specific. 
# In coming months we will investigate further whether there is value in additional cleanup stages, 
# most likely to be added into the pipeline after the ICA+FIX denoising. One approach that is simple but effective is “motion scrubbing”, #
# in which one identifies the timepoints that are “irreversibly” damaged by motion, and simply excises those from the timeseries analysis (Power et al., 2011).
#  We will evaluate this and other approaches, and where appropriate, make improvements to the temporal pre-processing pipeline.

if [ "$DO_ICA" -eq "1" ]; then
	printf "[$SUB] Performing MELODIC ICA ... "
	START=$(date -u +%s.%N)

	{

	DOIC=1
	if [ "$DOIC" -eq "1" ]; then	
		# Reorder columns in the motion estimates generated by AFNI to FSL's ordering
		touch ${OUTDIR}/motion_estimate_fsl.par
		while IFS=" " read -r roll pitch yaw ds dl dy
		do    
			rollRad=`echo ${roll} "*3.14159/180" | bc -l | awk '{printf "%f", $0}'`
			pitchRad=`echo ${pitch} "*3.14159/180" | bc -l | awk '{printf "%f", $0}'`
			yawRad=`echo ${yaw} "*3.14159/180" | bc -l | awk '{printf "%f", $0}'`

			echo $rollRad ' ' $pitchRad ' ' $yawRad ' ' $ds ' ' $dl ' ' $dy >> ${OUTDIR}/motion_estimate_fsl.par

		done < ${OUTDIR}/motion_estimate.par


		# Create a temporary detrended 3d+time series to run melodic
		# For ICA-FIX, smoothing is not recommended [only during visualization]
		3dTproject -prefix ${OUTDIR}/tmp_melodic.nii -polort 1 -stopband 0 0.009 -mask ${OUTDIR}/nat_mask.nii  -TR ${TR} -input ${OUTDIR}/proc_data_native.nii 

		if test -f "${OUTDIR}/func2anat0Warp.nii.gz"; then
			${ABIN}/antsApplyTransforms -v 1 \
				-i ${OUTDIR}/t1_frn.nii \
				--float \
				-r ${OUTDIR}/ref_func_bf.nii \
				-o ${OUTDIR}/t1_nat_bg.nii \
				-t [${OUTDIR}/func2anat.mat,1] \
				-t [${OUTDIR}/func2anat0InverseWarp.nii.gz,0] \
				-n Linear
		else
			${ABIN}/antsApplyTransforms -v 1 \
				-i ${OUTDIR}/t1_frn.nii \
				--float \
				-r ${OUTDIR}/ref_func_bf.nii \
				-o ${OUTDIR}/t1_nat_bg.nii \
				-t [${OUTDIR}/func2anat.mat,1] \
				-n Linear
		fi

		rm -f -R ${OUTDIR}/melodic.ic
		rm -f -R ${OUTDIR}/FIX

		# --dimest=lap is the default, but it seems empirically to overestimate components 
		# see also Varoquaux et al. 2010 NeuroImage
		melodic -i ${OUTDIR}/tmp_melodic.nii -o ${OUTDIR}/melodic.ic --tr=${TR} --mmthresh=0.5 --nobet --dimest=bic --approach=tica \
			 --mask=${OUTDIR}/nat_mask.nii --Ostats --report  --eps=0.0005 --bgimage=${OUTDIR}/t1_nat_bg.nii


		rm ${OUTDIR}/tmp_melodic.nii 

		# Creates a smoothed copy of the melodic components. Sometimes helps with visualization
		3dmerge -doall -1blur_fwhm 4 -prefix ${OUTDIR}/melodic.ic/melodic_IC_smooth.nii ${OUTDIR}/melodic.ic/melodic_IC.nii.gz


		# Creates a folder to be used with FIX
		# TODO avoid creating this if ICA-FIX procedure is not going to be used
		mkdir ${OUTDIR}/FIX
		mkdir ${OUTDIR}/FIX/mc
		mkdir ${OUTDIR}/FIX/reg
		mkdir ${OUTDIR}/FIX/filtered_func_data.ica

		cp ${OUTDIR}/proc_data_native.nii ${OUTDIR}/FIX/filtered_func_data.nii
		gzip ${OUTDIR}/FIX/filtered_func_data.nii
		cp ${OUTDIR}/motion_estimate_fsl.par ${OUTDIR}/FIX/mc/prefiltered_func_data_mcf.par
		cp -R ${OUTDIR}/melodic.ic/* ${OUTDIR}/FIX/filtered_func_data.ica
		cp ${OUTDIR}/nat_mask.nii ${OUTDIR}/FIX/mask.nii
		gzip ${OUTDIR}/FIX/mask.nii
		fslmaths ${OUTDIR}/proc_data_native.nii -Tmean ${OUTDIR}/FIX/mean_func
		cp ${OUTDIR}/FIX/mean_func.nii.gz ${OUTDIR}/FIX/reg/example_func.nii.gz
		cp ${OUTDIR}/t1_frn.nii ${OUTDIR}/FIX/reg/highres.nii
		gzip ${OUTDIR}/FIX/reg/highres.nii

		#~/Documents/itk/c3d/bin/c3d_affine_tool -ref ${OUTDIR}/ref_func_bf.nii -src ${OUTDIR}/t1_frn.nii ${OUTDIR}/func2anat.mat -inv -oitk ${OUTDIR}/FIX/highres2example_func.mat
		~/Documents/itk/c3d/bin/c3d_affine_tool ${OUTDIR}/afunc2anat.mat -inv -o ${OUTDIR}/iafunc2anat.mat
		cp ${OUTDIR}/iafunc2anat.mat ${OUTDIR}/FIX/reg/highres2example_func.mat

		

		CURDIR=`pwd`
		cd ${FIXBIN}
		source ${FIXBIN}/fix -f ${OUTDIR}/FIX
		cd ${CURDIR}

		rm -r -f  ${OUTDIR}/melodic.ic

	fi


	} &> ${OUTDIR}/04_ICA.log

	END=$(date -u +%s.%N)
	DIFF=`echo "( $END - $START )" | bc`
	printf "DONE [%.1f s]\n" $DIFF

fi

exit 1

# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------





# Apply normalisation WARP
if [ "$DO_NORM" -eq "1" ]; then

	printf "[$SUB] Warping functional images to MNI space ... "
	START=$(date -u +%s.%N)

	{

		if [ "${REG_MODEL}" == "FIX" ]; then

			rm -f ${OUTDIR}/proc_data_native_fix.nii
			rm -f ${OUTDIR}/proc_data_native_fixp.nii 

			CURDIR=`pwd`
			cd ${FIXBIN}
			source ${FIXBIN}/fix -c ${OUTDIR}/FIX ${FIX_CLASSIFIER} ${FIX_THR}
			# -A
			source ${FIXBIN}/fix -a ${OUTDIR}/FIX/fix4melview_${FIX_CL_LABEL}_thr${FIX_THR}.txt  
			cd ${CURDIR}

			mv ${OUTDIR}/FIX/filtered_func_data_clean.nii.gz ${OUTDIR}/proc_data_native_fix.nii.gz
			gunzip -f ${OUTDIR}/proc_data_native_fix.nii.gz


			if [ "${DO_PEST}" -eq "2" ]; then
				matlab  "-nodesktop -nosplash " <<<"apply_phycaa_rsn(${TR}, '${OUTDIR}', 'proc_data_native_fix', '${OUTDIR}/nat_mask.nii', '${OUTDIR}/csf_mask_nat.nii' ); exit;"
			fi
		fi



		norm_func()
		{

			TN=$1

			OUTDIR=$2
			ABIN=$3
			USBPATH=$4
			REF=${5}
			SMOOTH=${6}
			CORRECT_ANAT_FMAP_NATIVE=${7}

			if (( $TN < 10 )); then
				FI=00$TN
			elif (( $TN < 100 )); then
				FI=0$TN
			else
				FI=$TN	
			fi



			if [ "$SMOOTH" -eq "0" ]; then
				${ABIN}/DenoiseImage -s 1 -x ${OUTDIR}/nat_mask.nii -n Rician -v 1 -i ${OUTDIR}/tmp/proc_data.${FI}.nii  -o [ ${OUTDIR}/tmp/proc_data_s.${FI}.nii ]
				mv ${OUTDIR}/tmp/proc_data_s.${FI}.nii ${OUTDIR}/tmp/proc_data.${FI}.nii
			elif [ "$SMOOTH" -gt "0" ]; then
				# Smoothing in FSL is defined by sigma [see https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=fsl;54608a1b.1111]
				#You can use:
				#fslmaths original.nii -kernel gauss 2.1233226 -fmean smoothed.nii

				#The gaussian kernel takes its argument as sigma in mm instead of the FWHM, but 
				#you can see how to convert between these two values here:
				#http://mathworld.wolfram.com/GaussianFunction.html

				#Basically, divide your FWHM (in mm) by 2.3548 to get sigma.
				SIGMA=`echo "${SMOOTH}/2.3548" | bc -l`
				fslmaths ${OUTDIR}/tmp/proc_data.${FI}.nii -kernel gauss ${SIGMA} -fmean ${OUTDIR}/tmp/proc_data_s.${FI}.nii

				# NOTE: It is better to set FSL to save nii files by default (in terms of speed), though this might create problems for FEAT/FLAME (and certainly for FIX)
				# I have not tested FSL's pipeline after changing /usr/local/fsl/etc/fslconf/fsl.sh and /usr/local/fsl/etc/fslconf/fsl.csh [OUTPUT type]
				# Likely, this can be done at the start of this script, which is a TODO
				gunzip -f ${OUTDIR}/tmp/proc_data_s.${FI}.nii.gz
				mv ${OUTDIR}/tmp/proc_data_s.${FI}.nii ${OUTDIR}/tmp/proc_data.${FI}.nii
			fi


 
			if [ "${CORRECT_ANAT_FMAP_NATIVE}" -eq "0" ] && test -f "${OUTDIR}/func2anat0Warp.nii.gz"; then
				${ABIN}/antsApplyTransforms -v 0 \
					-i ${OUTDIR}/tmp/proc_data.$FI.nii \
					--float \
					-r ${REF} \
					-o ${OUTDIR}/tmp/proc_data_MNI_$FI.nii \
					-t [${OUTDIR}/anat2group1Warp.nii.gz,0] \
					-t [${OUTDIR}/anat2group0GenericAffine.mat,0] \
					-t [${OUTDIR}/func2anat0Warp.nii.gz,0] \
					-t [${OUTDIR}/func2anat.mat,0] \
					-n BSpline
			else
				${ABIN}/antsApplyTransforms -v 0 \
					-i ${OUTDIR}/tmp/proc_data.$FI.nii \
					--float \
					-r ${REF} \
					-o ${OUTDIR}/tmp/proc_data_MNI_$FI.nii \
					-t [${OUTDIR}/anat2group1Warp.nii.gz,0] \
					-t [${OUTDIR}/anat2group0GenericAffine.mat,0] \
					-t [${OUTDIR}/func2anat.mat,0] \
					-n BSpline
			fi

			# Remove negative values due to BSpline interpolation (mostly in the background of the images)
			3dmerge -prefix ${OUTDIR}/tmp/proc_data_MNI_thr_$FI.nii -1noneg ${OUTDIR}/tmp/proc_data_MNI_$FI.nii

		}
		export -f norm_func


		SMOOTH=-1
		if [ "$DO_NLM" -eq "1" ]; then
			print_debug 'Nonlocal Means filtering + Normalization'
			SMOOTH=0
			
		fi

	

		# In case, for some reason, this was not set before
		NVOLS=`fslnvols ${OUTDIR}/proc_data_native.nii`

		mkdir ${OUTDIR}/tmp
		3dTsplit4D -prefix ${OUTDIR}/tmp/proc_data.nii -keep_datum ${OUTDIR}/proc_data_native.nii
		parallel -j5 --line-buffer norm_func ::: $(seq 0 ${NVOLS}) ::: ${OUTDIR} ::: ${ABIN} ::: ${USBPATH} ::: ${MNI_REF_2mm} ::: $SMOOTH ::: ${CORRECT_ANAT_FMAP_NATIVE}
		rm ${OUTDIR}/proc_data_mni.nii
		3dTcat -prefix ${OUTDIR}/proc_data_mni.nii -tr ${TR} ${OUTDIR}/tmp/proc_data_MNI_thr_*.nii
		rm -r ${OUTDIR}/tmp

		# Perhaps only do for the desired level of removal (aggressive or not)
		if [ "${REG_MODEL}" == "FIX" ]; then


			mkdir ${OUTDIR}/tmp
			3dTsplit4D -prefix ${OUTDIR}/tmp/proc_data.nii -keep_datum ${OUTDIR}/proc_data_native_fix.nii
			parallel -j5 --line-buffer norm_func ::: $(seq 0 ${NVOLS}) ::: ${OUTDIR} ::: ${ABIN} ::: ${USBPATH} ::: ${MNI_REF_2mm} ::: $SMOOTH ::: ${CORRECT_ANAT_FMAP_NATIVE}
			rm ${OUTDIR}/proc_data_mni_fix.nii
			3dTcat -prefix ${OUTDIR}/proc_data_mni_fix.nii -tr ${TR} ${OUTDIR}/tmp/proc_data_MNI_thr_*.nii
			rm -r ${OUTDIR}/tmp

			if [ "${DO_PEST}" -eq "2" ]; then
				mkdir ${OUTDIR}/tmp
				3dTsplit4D -prefix ${OUTDIR}/tmp/proc_data.nii -keep_datum ${OUTDIR}/proc_data_native_fixp.nii
				parallel -j5 --line-buffer norm_func ::: $(seq 0 ${NVOLS}) ::: ${OUTDIR} ::: ${ABIN} ::: ${USBPATH} ::: ${MNI_REF_2mm} ::: $SMOOTH ::: ${CORRECT_ANAT_FMAP_NATIVE}
				rm ${OUTDIR}/proc_data_mni_fixp.nii
				3dTcat -prefix ${OUTDIR}/proc_data_mni_fixp.nii -tr ${TR} ${OUTDIR}/tmp/proc_data_MNI_thr_*.nii
				rm -r ${OUTDIR}/tmp
			fi

		fi
		3dAutomask -prefix ${OUTDIR}/mni_mask.nii ${OUTDIR}/proc_data_mni.nii



} &> ${OUTDIR}/05_norm2mni.log

	END=$(date -u +%s.%N)
	DIFF=`echo "( $END - $START )" | bc`
	printf "DONE [%.1f s]\n" $DIFF
fi





# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------

if [ "$DO_QA" -eq "1" ]; then

	printf "[$SUB] Generating QC plots ... "
	START=$(date -u +%s.%N)

{


#	python3.6 /home/fsluser/Documents/rs_proc/QC_funcs/QC_motion_estimate.py -out ${OUTDIR}/ -in ${OUTDIR}  -mpe motion_estimate.par \
#				-prog AFNI -outf Motion_Estimate -rngT 1.5 -rngR 1 -dpi 300 -tr ${TR}
#
#
#	python3.6 /home/fsluser/Documents/rs_proc/QC_funcs/QC_check_reg.py -out ${OUTDIR} -in ${OUTDIR}   \
#				-im1 ${OUTDIR}/t1_frn.nii -im2 ${OUTDIR}/afunc2anat.nii -outf 'func2anat' \
#				-dpi 300
#
#	fslmaths ${OUTDIR}/proc_data_mni.nii -Tmean ${OUTDIR}/mean_func_mni.nii
#	gunzip -f ${OUTDIR}/mean_func_mni.nii.gz
#	python3.6 /home/fsluser/Documents/rs_proc/QC_funcs/QC_check_reg.py -out ${OUTDIR} -in ${OUTDIR}   \
#				-im1 ${MNI_REF_2mm} -im2 ${OUTDIR}/mean_func_mni.nii -outf 'func2mni' \
#				-dpi 300


	fslmaths ${OUTDIR}/proc_data_mni.nii -Tmean ${OUTDIR}/mni_mean.nii
	gunzip -f ${OUTDIR}/mni_mean.nii.gz
	model_qa()
	{
		OUTDIR=${1}
		MNI_REF_2mm=${2}
		TR=${3}
		REG_MODEL=${4}
		mkdir ${OUTDIR}/QA_${REG_MODEL}
		QA_DIR=${OUTDIR}/QA_${REG_MODEL}

		rm -f ${OUTDIR}/proc_data_mni_${REG_MODEL}.nii
		rm -f ${OUTDIR}/proc_data_mni_${REG_MODEL}.nii.gz
		rm -f ${QA_DIR}/*

		
		if [ "${REG_MODEL}" == "NONE" ]; then
			3dTproject -prefix ${OUTDIR}/proc_data_mni_${REG_MODEL}.nii -polort 1 -mask ${OUTDIR}/mni_mask.nii -TR ${TR}  \
					 -stopband 0 0.01 \
					 -input ${OUTDIR}/proc_data_mni.nii	
		fi

		
		if [ "${REG_MODEL}" == "SFIX" ]; then

			3dTproject -prefix ${OUTDIR}/proc_data_mni_${REG_MODEL}.nii -polort 1 -mask ${OUTDIR}/mni_mask.nii -TR ${TR}  -cenmode NTRP \
					 -stopband 0 0.01  \
					 -ort ${OUTDIR}/motion_regressors_12.txt \
					 -ort ${OUTDIR}/motion_regressors_24.txt \
					 -censor ${OUTDIR}/temporal_mask_fd.txt \
					 -input ${OUTDIR}/proc_data_mni_fix.nii	
		fi

		if [ "${REG_MODEL}" == "SFIX_P" ]; then

			3dTproject -prefix ${OUTDIR}/proc_data_mni_${REG_MODEL}.nii -polort 1 -mask ${OUTDIR}/mni_mask.nii -TR ${TR}  -cenmode NTRP \
					 -stopband 0 0.01  \
					 -ort ${OUTDIR}/motion_regressors_12.txt \
					 -ort ${OUTDIR}/motion_regressors_24.txt \
					 -censor ${OUTDIR}/temporal_mask_fd.txt \
					 -input ${OUTDIR}/proc_data_mni_fixp.nii	
		fi

		if [ "${REG_MODEL}" == "SRP24WM1CSF1" ]; then
			3dTproject -prefix ${OUTDIR}/proc_data_mni_${REG_MODEL}.nii -polort 1  -mask ${OUTDIR}/mni_mask.nii -TR ${TR} -cenmode NTRP \
					 -stopband 0 0.01  \
					 -ort ${OUTDIR}/motion_regressors_12.txt \
					 -ort ${OUTDIR}/motion_regressors_24.txt \
					 -ort ${OUTDIR}/csf_sig.txt \
					 -ort ${OUTDIR}/wm_sig.txt \
					 -censor ${OUTDIR}/temporal_mask_fd.txt \
					 -input ${OUTDIR}/proc_data_mni.nii	

		fi

		
		if [ "${REG_MODEL}" == "SRP9" ]; then
			3dTproject -prefix ${OUTDIR}/proc_data_mni_${REG_MODEL}.nii -polort 1 -mask ${OUTDIR}/mni_mask.nii -TR ${TR} -cenmode NTRP \
					 -stopband 0 0.01   \
					 -ort ${OUTDIR}/motion_estimate.par \
					 -ort ${OUTDIR}/csf_sig.txt \
					 -ort ${OUTDIR}/wm_sig.txt \
					 -ort ${OUTDIR}/global_sig.txt \
					 -censor ${OUTDIR}/temporal_mask_fd.txt \
					 -input ${OUTDIR}/proc_data_mni.nii	

		fi

	

		if [ "${REG_MODEL}" == "SRP24CC" ]; then
			3dTproject -prefix ${OUTDIR}/proc_data_mni_${REG_MODEL}.nii -polort 1 -mask ${OUTDIR}/mni_mask.nii -TR ${TR} -cenmode NTRP \
					 -stopband 0 0.01  \
					 -ort ${OUTDIR}/motion_regressors_12.txt \
					 -ort ${OUTDIR}/motion_regressors_24.txt \
					 -ort ${OUTDIR}/acompcor.txt \
					 -censor ${OUTDIR}/temporal_mask_fd.txt \
					 -input ${OUTDIR}/proc_data_mni.nii	

		fi

		# TODO Only perform this step if intensity has not been glboally normalized previously
		#fslmaths ${OUTDIR}/proc_data_mni_${REG_MODEL}.nii -add ${OUTDIR}/mni_mean.nii -ing 1000 ${OUTDIR}/proc_data_mni_${REG_MODEL}_i.nii
		#gunzip -f ${OUTDIR}/proc_data_mni_${REG_MODEL}_i.nii.gz
		#mv ${OUTDIR}/proc_data_mni_${REG_MODEL}_i.nii ${OUTDIR}/proc_data_mni_${REG_MODEL}.nii



		python3.6 /home/fsluser/Documents/rs_proc/QC_funcs/QC_grey_plot.py -out ${QA_DIR} -in ${OUTDIR}  -mpe motion_estimate.par \
					-fname proc_data_mni_${REG_MODEL}.nii -csf_name csf_mask_mni.nii -wm_name wm_mask_mni.nii -gm_name gm_mask_mni.nii \
					-norm zscore -range 1.0 \
					-prog AFNI -outf 02_Greyplot_${REG_MODEL}  -dpi 300 -tr ${TR}
#

#		python3.6 /home/fsluser/Documents/rs_proc/QC_funcs/QC_temporal_stats.py -out ${QA_DIR} -in ${OUTDIR}  -mpe motion_estimate.par \
#					-fname proc_data_mni_${REG_MODEL}.nii -bg ${MNI_REF_2mm} -type std \
#					-range 95% -plane axial -save_nii 1 \
#					-prog AFNI -outf 03_TemporalStd_z_${REG_MODEL}  -dpi 300 
#
#		python3.6 /home/fsluser/Documents/rs_proc/QC_funcs/QC_temporal_stats.py -out ${QA_DIR} -in ${OUTDIR}  -mpe motion_estimate.par \
#					-fname proc_data_mni_${REG_MODEL}.nii -bg ${MNI_REF_2mm} -type std \
#					-range 95% -plane coronal \
#					-prog AFNI -outf 03_TemporalStd_y_${REG_MODEL}  -dpi 300 
#
#		python3.6 /home/fsluser/Documents/rs_proc/QC_funcs/QC_temporal_stats.py -out ${QA_DIR} -in ${OUTDIR}  -mpe motion_estimate.par \
#					-fname proc_data_mni_${REG_MODEL}.nii -bg ${MNI_REF_2mm} -type std \
#					-range 95% -plane sagital \
#					-prog AFNI -outf 03_TemporalStd_x_${REG_MODEL}  -dpi 300 
#

#		python3.6 /home/fsluser/Documents/rs_proc/QC_funcs/QC_temporal_stats.py -out ${QA_DIR} -in ${OUTDIR}  -mpe motion_estimate.par \
#					-fname proc_data_mni_${REG_MODEL}.nii -bg ${MNI_REF_2mm} -type GlobalCorr \
#					-range 95% -plane coronal -thr 0.3 -smooth 0 \
#					-prog AFNI -outf 04_GlobalCorr_y_${REG_MODEL}  -dpi 300 
#
		python3.6 /home/fsluser/Documents/rs_proc/QC_funcs/QC_temporal_stats.py -out ${QA_DIR} -in ${OUTDIR}  -mpe motion_estimate.par \
					-fname proc_data_mni_${REG_MODEL}.nii -bg ${MNI_REF_2mm} -type GlobalCorr \
					-range 95% -plane axial -thr 0.3 -smooth 0 -save_nii 1 \
					-prog AFNI -outf 04_GlobalCorr_z_${REG_MODEL}  -dpi 300 
#

#		python3.6 /home/fsluser/Documents/rs_proc/QC_funcs/QC_temporal_stats.py -out ${QA_DIR} -in ${OUTDIR}  -mpe motion_estimate.par \
#					-fname proc_data_mni_${REG_MODEL}.nii -bg ${MNI_REF_2mm} -type GlobalCorr \
#					-range 95% -plane sagital -thr 0.3 -smooth 0 \
#					-prog AFNI -outf 04_GlobalCorr_x_${REG_MODEL}  -dpi 300 
#
#
		python3.6 /home/fsluser/Documents/rs_proc/QC_funcs/QC_temporal_stats.py -out ${QA_DIR} -in ${OUTDIR}  -mpe motion_estimate.par \
					-fname proc_data_mni_${REG_MODEL}.nii -bg ${MNI_REF_2mm} -type MotionCorr \
					-range 95% -plane axial -thr 0.3 -smooth 0 -save_nii 1 \
					-prog AFNI -outf 04_MotionCorr_z_${REG_MODEL}  -dpi 300 
#
#		python3.6 /home/fsluser/Documents/rs_proc/QC_funcs/QC_temporal_stats.py -out ${QA_DIR} -in ${OUTDIR}  -mpe motion_estimate.par \
#					-fname proc_data_mni_${REG_MODEL}.nii -bg ${MNI_REF_2mm} -type MotionCorr \
#					-range 95% -plane coronal -thr 0.3 -smooth 0 \
#					-prog AFNI -outf 04_MotionCorr_y_${REG_MODEL}  -dpi 300 
#
#		python3.6 /home/fsluser/Documents/rs_proc/QC_funcs/QC_temporal_stats.py -out ${QA_DIR} -in ${OUTDIR}  -mpe motion_estimate.par \
#					-fname proc_data_mni_${REG_MODEL}.nii -bg ${MNI_REF_2mm} -type MotionCorr \
#					-range 95% -plane sagital -thr 0.3 -smooth 0 \
#					-prog AFNI -outf 04_MotionCorr_x_${REG_MODEL}  -dpi 300 


		

#		python3.6 /home/fsluser/Documents/rs_proc/QC_funcs/QC_FC.py -out ${QA_DIR} -in ${OUTDIR}   -tr ${TR} \
#				-fname proc_data_mni_${REG_MODEL}.nii -outf 05_FC_${REG_MODEL}  \
#				-low 0.1 -high 0.01 -dpi 300 

		python3.6 /home/fsluser/Documents/rs_proc/QC_funcs/QC_FC.py -out ${QA_DIR} -in ${OUTDIR}   -tr ${TR} \
				-fname proc_data_mni_${REG_MODEL}.nii -outf 05_FC_${REG_MODEL}  \
				 -dpi 300 

		#rm -f ${OUTDIR}/proc_data_mni_${REG_MODEL}.nii
		#rm -f ${OUTDIR}/proc_data_mni_${REG_MODEL}.nii.gz
	}
	export -f model_qa

	parallel -j3 --line-buffer model_qa ::: ${OUTDIR} ::: ${MNI_REF_2mm} ::: ${TR} :::  NONE SFIX SFIX_P SRP24WM1CSF1 SRP9 SRP24CC #NONE SFIXM SRP24WM1CSF1 SFIX_NONAGG SRP9 SRP24CC #NONE FIX_NONAGG FIX_AGG SFIX_NONAGG SFIX_AGG  RP24WM1CSF1 RP9 RP24CC SRP24WM1CSF1 SRP9 SRP24CC



} &> ${OUTDIR}/06_qa.log

	END=$(date -u +%s.%N)
	DIFF=`echo "( $END - $START )" | bc`
	printf "DONE [%.1f s]\n" $DIFF
fi

#cp -R /mnt/hgfs/CRUNCHOUT/RS${SUB}/* ${OUTDIR}/
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------

exit 1


if [ "$DO_CLEAN" -eq "1" ]; then

	#=====================================



	{

		rm -f ${OUTDIR}/func_data_*.nii*

	} &> /dev/null



fi

END_ALL=$(date -u +%s.%N)
DIFF=`echo "( $END_ALL - $START_ALL )" | bc`
printf "DONE in %.1f seconds\n" $DIFF 

exit 1

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






