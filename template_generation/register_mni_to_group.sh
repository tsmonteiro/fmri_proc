#!/bin/bash

# Example call
# ./import_t1s.sh SUB001 import_sub_t1.sh /path/to/all_anats/
ABIN=/home/luna.kuleuven.be/u0101486/ANTs/bin/ # PATH where ANTs is installed
WDIR=/home/luna.kuleuven.be/u0101486/workspace/fmri_proc/
ABIN=/home/luna.kuleuven.be/u0101486/ANTs/bin/ # PATH where ANTs is installed

#GROUPDIR='/home/luna.kuleuven.be/u0101486/workspace/data/RepImpact/Template/'
#GROUPDIR='/home/luna.kuleuven.be/u0101486/workspace/data/ConnectEx/Template/'
#GROUPDIR='/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/Templates/RSN/Mantini'
GROUPDIR='/home/luna.kuleuven.be/u0101486/workspace/data/RSPET/Template/'

# ONCE The group average is generated, run the code below
rm -f ${GROUPDIR}/template_mni.nii

#python ${WDIR}/extract_brain.py -i "${GROUPDIR}/Group_T1_Avg_Brain.nii" -o "${GROUPDIR}/brain_mask.nii"
#3dcalc -a ${GROUPDIR}/Group_T1_Avg_Brain.nii -b ${GROUPDIR}/brain_mask.nii -expr "a*b" -prefix ${GROUPDIR}/Group_T1_Avg_Brain_b.nii

# Register MNI to GROUP
cp ${WDIR}/atlas/CAT12/lg400_cobra_IXI.nii ${GROUPDIR}/lg400_cobra_mni.nii
cp ${WDIR}/atlas/CAT12/yeo17cobra_IXI.nii ${GROUPDIR}/yeo17cobra_mni.nii
cp ${WDIR}/atlas/CAT12/FSL_MNI152_IXI_Brain.nii ${GROUPDIR}/template_mni.nii

# For ConnectEx, copy the Mantini 2013 template and networks
#cp ${WDIR}/atlas/Mantini2013/ch2.nii ${GROUPDIR}/template_mni.nii

#NII_IN=(cingulo_opercular default_mode dorsal_attention dorsal_somatomotor \
#      ventral_attention early_auditory language lateral_prefrontal \
#      left_frontoparietal medial_prefrontal right_frontoparietal \
#      ventral_somatomotor visual_parafoveal visual_peripheral)

#for NII in "${NII_IN[@]}";
#do
#  cp ${WDIR}/atlas/Mantini2013/RSNs_fMRI/${NII}.nii ${GROUPDIR}/${NII}_mni.nii
#done


REF=${GROUPDIR}/Group_T1_Avg_Brain.nii
SRC=${GROUPDIR}/template_mni.nii
#REF=/home/luna.kuleuven.be/u0101486/workspace/fmri_proc/atlas/CAT12/FSL_MNI152_FreeSurferConformed_1mm.nii
#SRC=/home/luna.kuleuven.be/u0101486/workspace/fmri_proc/atlas/CAT12/Template_T1_IXI555_MNI152_GS.nii
${ABIN}/antsRegistration -d 3 -r [$REF,$SRC,0] -v 1 \
            -m MI[$REF,$SRC,1,32] -t translation[0.1] -c [500,5.e-7,20] \
            -s 3vox -f 3 -l 1 -n BSpline \
            -m MI[$REF,$SRC,1,32,Regular,0.25] -t rigid[0.1] -c [500,5.e-7,20] \
            -s 3vox -f 3 -l 1 -n BSpline \
            -m MI[$REF,$SRC,1,32,Regular,0.25] -t affine[0.1] -c [500x100x10,5.e-7,10] \
            -s 2x1x0vox -f 3x2x1 -l 1 -n BSpline \
            -m CC[$REF,$SRC,1,3] -t SyN[0.2,3] -c [30x10,1.e-7,10] \
            -s 1x0vox -f 2x1 -l 1 -n BSpline \
            -o [${GROUPDIR}/mni2group,${GROUPDIR}/mni2group.nii]
            #-o [/home/luna.kuleuven.be/u0101486/workspace/fmri_proc/atlas/CAT12/IXI_2_FSL,/home/luna.kuleuven.be/u0101486/workspace/fmri_proc/atlas/CAT12/IXI_2_FSL.nii]
#-o [${GROUPDIR}/mni2group,${GROUPDIR}/mni2group.nii]


for NII in "${NII_IN[@]}";
do
  cp ${WDIR}/atlas/Mantini2013/RSNs_fMRI/${NII}.nii ${GROUPDIR}/${NII}_mni.nii

  ${ABIN}/antsApplyTransforms \
        -i ${GROUPDIR}/${NII}_mni.nii \
        --float \
        -r ${GROUPDIR}/Group_T1_Avg_Brain.nii \
        -o ${GROUPDIR}/${NII}_group.nii \
        -t [${GROUPDIR}/mni2group0GenericAffine.mat,0] \
        -t [${GROUPDIR}/mni2group1Warp.nii.gz,0] \
        -n BSpline[5]

  fslmaths ${GROUPDIR}/${NII}_group.nii -thr 0 ${GROUPDIR}/${NII}_group_corr.nii
  gunzip -f ${GROUPDIR}/${NII}_group_corr.nii.gz
  mv ${GROUPDIR}/${NII}_group_corr.nii ${GROUPDIR}/${NII}_group.nii
done





#exit

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
