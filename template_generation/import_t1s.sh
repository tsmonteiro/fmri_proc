#!/bin/bash

# Example call
# ./import_t1s.sh SUB001 import_sub_t1.sh /path/to/all_anats/

SUB=$1
IMPORTFILE=$2 #/home/fsluser/Documents/rs_proc/conf_files/rep_impact_belgium.conf
OUTDIR=$3
ABIN=/media/thiago/EXTRALINUX/ANTs/bin/ # PATH where ANTs is installed

source ${IMPORTFILE} ${SUB} ${OUTDIR}

if [ ! -f "${OUTDIR}/t1_${SUB}.nii" ]; then
  # Failed to import T1
	exit
fi

# Correct the orientation of T1 to somethign close to standard orientation
matlab "-nodesktop -nosplash " <<<"cd './matlab'; p_reorient('${OUTDIR}/t1_${SUB}.nii'); exit;"

fslmaths ${OUTDIR}/t1_${SUB}.nii -thr 10 -nan ${OUTDIR}/t1_${SUB}_t.nii
gunzip -f ${OUTDIR}/t1_${SUB}_t.nii.gz

mv ${OUTDIR}/t1_${SUB}_t.nii ${OUTDIR}/t1_${SUB}.nii



## The following 2 lines should be only executed for subjects whose T1 are strongly tilted forward (-) or backward
## In those instances SPm segmentation will fail
## !!IMPORTANT!! The angle MUST BE INDIVIDUALLY ADJUSTED!!
#3drotate -prefix ${OUTDIR}/t1_${SUB}_rr.nii -rotate -10R 0 0  ${OUTDIR}/t1_${SUB}.nii
#mv ${OUTDIR}/t1_${SUB}_rr.nii ${OUTDIR}/t1_${SUB}.nii

3dresample -prefix ${OUTDIR}/t1_${SUB}_r.nii -orient ras -input ${OUTDIR}/t1_${SUB}.nii
3drefit -deoblique ${OUTDIR}/t1_${SUB}_r.nii

#3dAutomask -prefix ${OUTDIR}/brain_mask_${SUB}.nii -SI 130 ${OUTDIR}/t1_${SUB}_r.nii

${ABIN}/DenoiseImage -i ${OUTDIR}/t1_${SUB}_r.nii -s 2 -r 2 -p 1 -n Rician -o [ ${OUTDIR}/t1_${SUB}_rf.nii ]
fslmaths ${OUTDIR}/t1_${SUB}_rf.nii -thr 10 -nan ${OUTDIR}/t1_${SUB}_rft.nii
gunzip -f ${OUTDIR}/t1_${SUB}_rft.nii.gz
mv ${OUTDIR}/t1_${SUB}_rft.nii ${OUTDIR}/t1_${SUB}_rf.nii


rm -f ${OUTDIR}/t1_${SUB}_r.nii
rm -f ${OUTDIR}/t1_${SUB}.nii
matlab "-nodesktop -nosplash " <<<"cd template_generation; spm_dartel_segment('${OUTDIR}/t1_${SUB}_rf.nii'); exit;"


rm -f ${OUTDIR}/c3t1_${SUB}_rf.nii
rm -f ${OUTDIR}/c4t1_${SUB}_rf.nii
rm -f ${OUTDIR}/c5t1_${SUB}_rf.nii


# ONCE The group average is generated, run the code below

# Register MNI to GROUP
#cp ${WDIR}/atlas/CAT12/lg400_cobra_IXI.nii ${GROUPDIR}/lg400_cobra_mni.nii
#cp ${WDIR}/atlas/CAT12/yeo17cobra_IXI.nii ${GROUPDIR}/yeo17cobra_mni.nii
#cp ${WDIR}/atlas/CAT12/FSL_MNI152_IXI_Brain.nii ${GROUPDIR}/template_mni.nii
#
#REF=${GROUPDIR}/Group_T1_Avg_Brain.nii
#SRC=${GROUPDIR}/template_mni.nii
#${ABIN}/antsRegistration -d 3 -r [$REF,$SRC,0] -v 1 \
#            -m MI[$REF,$SRC,1,32] -t translation[0.1] -c [500,5.e-7,20] \
#            -s 3vox -f 3 -l 1 -n BSpline \
#            -m MI[$REF,$SRC,1,32,Regular,0.25] -t rigid[0.1] -c [500,5.e-7,20] \
#            -s 3vox -f 3 -l 1 -n BSpline \
#            -m MI[$REF,$SRC,1,32,Regular,0.25] -t affine[0.1] -c [500x100x10,5.e-7,10] \
#            -s 2x1x0vox -f 3x2x1 -l 1 -n BSpline \
#            -m CC[$REF,$SRC,1,3] -t SyN[0.2,3] -c [20x10,1.e-7,10] \
#            -s 1x0vox -f 2x1 -l 1 -n BSpline \
#            -o [${GROUPDIR}/mni2group,${GROUPDIR}/mni2group.nii]
#
#
#${ABIN}/antsApplyTransforms \
#      -i ${GROUPDIR}/lg400_cobra_mni.nii \
#      --float \
#      -r ${GROUPDIR}/Group_T1_Avg_Brain.nii \
#      -o ${GROUPDIR}/lg400_cobra_group.nii \
#      -t [${GROUPDIR}/mni2group0GenericAffine.mat,0] \
#      -t [${GROUPDIR}/mni2group1Warp.nii.gz,0] \
#      -n NearestNeighbor
#
#${ABIN}/antsApplyTransforms \
#      -i ${GROUPDIR}/yeo17cobra_mni.nii \
#      --float \
#      -r ${GROUPDIR}/Group_T1_Avg_Brain.nii \
#      -o ${GROUPDIR}/yeo17cobra_group.nii \
#      -t [${GROUPDIR}/mni2group0GenericAffine.mat,0] \
#      -t [${GROUPDIR}/mni2group1Warp.nii.gz,0] \
#      -n NearestNeighbor
