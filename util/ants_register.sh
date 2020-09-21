#!/bin/bash

# Example call
# ./import_t1s.sh SUB001 import_sub_t1.sh /path/to/all_anats/
ABIN=/home/luna.kuleuven.be/u0101486/ANTs/bin/ # PATH where ANTs is installed
WDIR=/home/luna.kuleuven.be/u0101486/workspace/fmri_proc/
ABIN=/home/luna.kuleuven.be/u0101486/ANTs/bin/ # PATH where ANTs is installed



BASEDIR='/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/iCAPs_results/dFC_Full3_Alpha_5_95_Fraction_0DOT05/K_15_Dist_cosine_Folds_17/Match_to_CC_2019/'
# ONCE The group average is generated, run the code below


#3dTsplit4D -prefix ${BASEDIR}/ICAP_3d/ICAP_Map.nii ${BASEDIR}/iCAPs_z.nii
#exit

REF=${BASEDIR}/Group_T1_Avg_Brain.nii
SRC=${BASEDIR}/ch2.nii
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
            -o [${BASEDIR}/mni2group,${BASEDIR}/mni2group.nii]



NII_IN=(00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 \
  21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 \
  43 44 45 46 47 48 49)

NII_IN=(cingulo_opercular default_mode dorsal_attention dorsal_somatomotor \
      ventral_attention early_auditory language lateral_prefrontal \
      left_frontoparietal medial_prefrontal right_frontoparietal \
      ventral_somatomotor visual_parafoveal visual_peripheral)

#NII_IN=(00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 )
for NII in "${NII_IN[@]}";
do

        #-t [${BASEDIR}/mni2group0GenericAffine.mat,0] \
        #-t [${BASEDIR}/mni2group1Warp.nii.gz,0] \
  ${ABIN}/antsApplyTransforms \
        -i ${BASEDIR}/mantini/${NII}.nii \
        --float \
        -r ${REF} \
         -t [${BASEDIR}/mni2group0GenericAffine.mat,0] \
         -t [${BASEDIR}/mni2group1Warp.nii.gz,0] \
        -o ${BASEDIR}/mantini/${NII}.Group.nii \
        -n BSpline[5]
done





exit
