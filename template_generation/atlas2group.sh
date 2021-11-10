#!/bin/bash

# Example call
# ./import_t1s.sh SUB001 import_sub_t1.sh /path/to/all_anats/


ABIN=/media/thiago/EXTRALINUX/ANTs/abin/bin/ # PATH where ANTs is installed
WDIR=/media/thiago/EXTRALINUX/fmri_proc/
GROUPDIR=/media/thiago/EXTRALINUX/fmri_proc/DATA/Template
# Register MNI to GROUP
cp ${WDIR}/atlas/CAT12/lg400_cobra_IXI.nii ${GROUPDIR}/lg400_cobra_mni.nii
cp ${WDIR}/atlas/CAT12/yeo17cobra_IXI.nii ${GROUPDIR}/yeo17cobra_mni.nii
cp ${WDIR}/atlas/CAT12/FSL_MNI152_IXI_Brain.nii ${GROUPDIR}/template_mni.nii

REF=${GROUPDIR}/Group_T1_Avg_Brain.nii
SRC=${GROUPDIR}/template_mni.nii
${ABIN}/antsRegistration -d 3 -r [$REF,$SRC,0] -v 1 \
            -m MI[$REF,$SRC,1,32] -t translation[0.1] -c [500,5.e-7,20] \
            -s 3vox -f 3 -l 1 -n BSpline \
            -m MI[$REF,$SRC,1,32,Regular,0.25] -t rigid[0.1] -c [500,5.e-7,20] \
            -s 3vox -f 3 -l 1 -n BSpline \
            -m MI[$REF,$SRC,1,32,Regular,0.25] -t affine[0.1] -c [500x100x10,5.e-7,10] \
            -s 2x1x0vox -f 3x2x1 -l 1 -n BSpline \
            -m CC[$REF,$SRC,1,3] -t SyN[0.2,3] -c [20x10,1.e-7,10] \
            -s 1x0vox -f 2x1 -l 1 -n BSpline \
            -o [${GROUPDIR}/mni2group,${GROUPDIR}/mni2group.nii]


${ABIN}/antsApplyTransforms \
      -i ${GROUPDIR}/lg400_cobra_mni.nii \
      --float \
      -r ${GROUPDIR}/Group_T1_Avg_Brain.nii \
      -o ${GROUPDIR}/lg400_cobra_group.nii \
      -t [${GROUPDIR}/mni2group0GenericAffine.mat,0] \
      -t [${GROUPDIR}/mni2group1Warp.nii.gz,0] \
      -n NearestNeighbor

${ABIN}/antsApplyTransforms \
      -i ${GROUPDIR}/yeo17cobra_mni.nii \
      --float \
      -r ${GROUPDIR}/Group_T1_Avg_Brain.nii \
      -o ${GROUPDIR}/yeo17cobra_group.nii \
      -t [${GROUPDIR}/mni2group0GenericAffine.mat,0] \
      -t [${GROUPDIR}/mni2group1Warp.nii.gz,0] \
      -n NearestNeighbor
