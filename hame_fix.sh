#!/bin/bash

# Corrects geometric distortion and regresses out nuisance covariates
./proc_rsdata.sh 054 /home/fsluser/Documents/rs_proc/conf_files/crunch.conf 

OUTDIR=$1



NVOLS=`fslnvols ${OUTDIR}/proc_data_native.nii`
native_fmap_corr()
{

	TN=$1

	OUTDIR=$2
	ABIN=$3
	USBPATH=$4
	REF=${OUTDIR}/ref_func_bfm.nii
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
		-o ${OUTDIR}/tmp/proc_data_MNI_$FI.nii \
		-t [${OUTDIR}/func2anat.mat,0] \
		-t [${OUTDIR}/func2anat0Warp.nii.gz,0] \
		-t [${OUTDIR}/func2anat.mat,1] \
		-n BSpline

	# Remove negative values due to BSpline interpolation (mostly in the background of the images)
	3dmerge -prefix ${OUTDIR}/tmp/proc_data_MNI_thr_$FI.nii -1noneg ${OUTDIR}/tmp/proc_data_MNI_$FI.nii

}
export -f native_fmap_corr



mkdir ${OUTDIR}/tmp
3dTsplit4D -prefix ${OUTDIR}/tmp/proc_data.nii -keep_datum ${OUTDIR}/proc_data_native.nii
parallel -j5 --line-buffer native_fmap_corr ::: $(seq 0 $NVOLS) ::: ${OUTDIR} ::: ${ABIN} ::: ${USBPATH}
rm ${OUTDIR}/proc_data_MNI2mm.nii

3dTcat -prefix ${OUTDIR}/proc_data_native_fmap.nii -tr ${TR} ${OUTDIR}/tmp/proc_data_MNI_thr_*.nii

rm -r ${OUTDIR}/tmp

# Enable this line later
#mv ${OUTDIR}/proc_data_native_fmap.nii ${OUTDIR}/proc_data_native_fmap.nii

3dTproject -prefix ${OUTDIR}/proc_data_native_clean.nii -polort 1 -mask ${OUTDIR}/nat_mask.nii -TR ${TR} -cenmode NTRP \
			 -stopband 0 0.01 \
			 -ort ${OUTDIR}/motion_regressors_12.txt \
			 -ort ${OUTDIR}/motion_regressors_24.txt \
			 -ort ${OUTDIR}/csf_sig.txt \
			 -ort ${OUTDIR}/wm_sig.txt \
			 -censor ${OUTDIR}/temporal_mask_fd.txt \
			 -input ${OUTDIR}/proc_data_native_fmap.nii	

rm ${OUTDIR}/proc_data_native_fmap.nii


