#!/bin/bash
# Clear terminal screen
#printf "\033c"


PIXDIM=`nifti_tool -quiet -disp_hdr -field pixdim -infiles ../data/Tests/func_data.nii`


#mkdir ../data/Tests/tmp/

IMBASE_O=../data/Tests//func_data.nii
IMBASE=../data/Tests//wfunc_data.nii
IMTARGET_O=../data/Tests/t1.nii
IMTARGET=../data/Tests/wt1.nii
IMTARGET_N=../data/Tests/nwt1.nii
IMTARGET_N2=../data/Tests/nafunc2anat.nii

rm -f ${IMTARGET_N}
rm -f ${IMTARGET_N2}
3dresample -prefix ${IMTARGET_N} -dxyz 2.375 2.375 2.75 -input ${IMTARGET} -master ${IMBASE}
3dresample -prefix ${IMTARGET_N2} -dxyz 2.375 2.375 2.75 -input ../data/Tests/afunc2anat.nii -master ${IMBASE}

#3dresample -prefix ${IMBASE} -orient ras -input ${IMBASE_O}
#3dresample -prefix ${IMTARGET} -orient ras -input ${IMTARGET_O}

#3dTsplit4D -prefix ../data/Tests/tmp/func_data.nii -keep_datum ${IMBASE}
IMBASE_3D=../data/Tests/tmp/func_data.000.nii
#rm -f ../data/Tests/afunc2anat.nii
#3dAllineate -input ${IMBASE_3D} -base ${IMTARGET_N} -source_automask -automask \
#    -prefix ../data/Tests/afunc2anat.nii -1Dmatrix_save ../data/Tests/init_alignment \
#    -1Dparam_save ../data/Tests/init_params -onepass -cost lpc+hel -warp shift_rotate

rm -f ../data/Tests/rfunc_data.nii
3dvolreg -heptic -prefix ../data/Tests/rfunc_data.nii -base ${IMTARGET_N} -rot_thresh 0.01 -delta 2 \
      -x_thresh 0.01 -zpad 10 -maxite 10 ${IMBASE}


rm -f ../data/Tests/r2func_data.nii
3dvolreg -heptic -prefix ../data/Tests/r2func_data.nii -base ${IMTARGET_N2} -rot_thresh 0.01 -delta 2 \
      -x_thresh 0.01 -zpad 10 -maxite 10 ${IMBASE}
