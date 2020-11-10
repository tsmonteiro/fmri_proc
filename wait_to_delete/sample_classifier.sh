#!/bin/bash
# Clear terminal screen
printf "\033c"

# Change as needed for each project
ABIN=/home/luna.kuleuven.be/u0101486/ANTs/bin/ # PATH where ANTs is installed
FIXBIN=/home/luna.kuleuven.be/u0101486/Software/fix/ # PATH where ICA_FIX is installed
WDIR='/home/luna.kuleuven.be/u0101486/workspace/fmri_proc/'
#BASEDIR=/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/
BASEDIR=/home/luna.kuleuven.be/u0101486/workspace/data/RepImpact/tmp/
#BASEDIR=/home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/tmp/
MNI_REF=/home/luna.kuleuven.be/u0101486/workspace/data/RepImpact/Template/Template_ 		# DARTEL Template

#CLASSIFIER=/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/fix_classifier_rs/classifier
CLASSIFIER=/home/luna.kuleuven.be/u0101486/workspace/data/RepImpact/fix_sample_n/
#CLASSIFIER=/home/luna.kuleuven.be/u0101486/workspace/data/RepImpact/fix_classifier_norway/classifier
#CLASSIFIER=/home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/fix_classifier_rs/classifier


# CRUNCH
#SUBS=(${BASEDIR}/RS087/FIX ${BASEDIR}/RS020/FIX ${BASEDIR}/RS086/FIX ${BASEDIR}/RS079/FIX \
#      ${BASEDIR}/RS075/FIX ${BASEDIR}/RS006/FIX ${BASEDIR}/RS055/FIX ${BASEDIR}/RS080/FIX \
#      ${BASEDIR}/RS093/FIX ${BASEDIR}/RS008/FIX ${BASEDIR}/RS028/FIX ${BASEDIR}/RS031/FIX \
#      ${BASEDIR}/RS037/FIX ${BASEDIR}/RS024/FIX ${BASEDIR}/RS022/FIX)

# RepImpact BELGIUM
SUBS=(B3_17 B2_30 B2_11 \
      B3_37 B3_41 B3_32 B1_07 \
      B3_58 B1_57 B2_48 B1_63 \
      B3_30 B2_60 B1_23 B1_24 \
      B1_09 B2_26 )

# RepImpact NORWAY
SUBS=(N1_06 N1_14 N1_21 N1_27 \
      N1_34 N1_52 N2_43 N2_05 \
      N2_26 N2_03 N2_40 N2_37 \
      N3_06 N3_15 N3_39 N3_02 \
      N3_18 )


# CAI China
#SUBS=(${BASEDIR}/sub101/FIX ${BASEDIR}/sub023/FIX ${BASEDIR}/sub008/FIX ${BASEDIR}/sub115/FIX \
#      ${BASEDIR}/sub122/FIX ${BASEDIR}/sub018/FIX ${BASEDIR}/sub015/FIX ${BASEDIR}/sub110/FIX \
#      ${BASEDIR}/sub109/FIX ${BASEDIR}/sub012/FIX ${BASEDIR}/sub004/FIX ${BASEDIR}/sub124/FIX)
NSUBS=0




for SUB in ${SUBS[@]}; do
  NSUBS=$((NSUBS+1))
  3dAutomask -prefix ${CLASSIFIER}/mask_${SUB}.nii ${CLASSIFIER}/${SUB}_group.nii
done

3dTcat -prefix ${CLASSIFIER}/all_masks.nii ${CLASSIFIER}/mask_*

exit

SI=0
#for SUB in ${SUBS[@]}; do
function norm_sub(){

      SUB=$1
      ABIN=/home/luna.kuleuven.be/u0101486/ANTs/bin/ # PATH where ANTs is installed
      FIXBIN=/home/luna.kuleuven.be/u0101486/Software/fix/ # PATH where ICA_FIX is installed
      WDIR='/home/luna.kuleuven.be/u0101486/workspace/fmri_proc/'
      #BASEDIR=/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/
      BASEDIR=/home/luna.kuleuven.be/u0101486/workspace/data/RepImpact/tmp/
      #BASEDIR=/home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/tmp/
      MNI_REF=/home/luna.kuleuven.be/u0101486/workspace/data/RepImpact/Template/Template_ 		# DARTEL Template
      ANATDIR=${BASEDIR}/${SUB}/anat

      #CLASSIFIER=/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/fix_classifier_rs/classifier
      CLASSIFIER=/home/luna.kuleuven.be/u0101486/workspace/data/RepImpact/fix_sample_n/
      cp ${BASEDIR}/${SUB}/proc_data_native_nlm.nii ${CLASSIFIER}/${SUB}.nii


      3dresample -prefix ${CLASSIFIER}/${SUB}_tmp.nii -orient ras -input ${CLASSIFIER}/${SUB}.nii
      3drefit -deoblique ${CLASSIFIER}/${SUB}_tmp.nii


      cp -f ${ANATDIR}/anat_proc.nii ${CLASSIFIER}/${SUB}_T1.nii
      cp -f ${ANATDIR}/anat_proc_brain.nii ${CLASSIFIER}/${SUB}_T1_brain.nii
      3dresample -prefix ${CLASSIFIER}/${SUB}_T1_tmp.nii -orient ras -input ${CLASSIFIER}/${SUB}_T1.nii
      3drefit -deoblique ${CLASSIFIER}/${SUB}_T1_tmp.nii

      3dresample -prefix ${CLASSIFIER}/${SUB}_T1_brain_tmp.nii -orient ras -input ${CLASSIFIER}/${SUB}_T1_brain.nii
      3drefit -deoblique ${CLASSIFIER}/${SUB}_T1_brain_tmp.nii

      mv ${CLASSIFIER}/${SUB}_tmp.nii ${CLASSIFIER}/${SUB}.nii
      mv ${CLASSIFIER}/${SUB}_T1_tmp.nii ${CLASSIFIER}/${SUB}_T1.nii
      mv ${CLASSIFIER}/${SUB}_T1_brain_tmp.nii ${CLASSIFIER}/${SUB}_T1_brain.nii

      SI=$((SI+1))


      NVOLS=`fslnvols ${CLASSIFIER}/${SUB}.nii`


      #FUNC
      fslmaths ${CLASSIFIER}/${SUB}.nii -Tmean -thrp 30 ${CLASSIFIER}/${SUB}_mean.nii
      gunzip -f ${CLASSIFIER}/${SUB}_mean.nii.gz
      ${ABIN}/DenoiseImage -i ${CLASSIFIER}/${SUB}_mean.nii -r 3x3x3 -n Gaussian -o [ ${CLASSIFIER}/${SUB}_mean_nds.nii ]



      copy_header ${CLASSIFIER}/${SUB}_mean.nii ${CLASSIFIER}/${SUB}_mean_nds.nii

      matlab "-nodesktop -nosplash " <<<"coreg_normalise('${CLASSIFIER}/${SUB}_mean_nds.nii', '${CLASSIFIER}/${SUB}.nii', ${NVOLS}, '${CLASSIFIER}/${SUB}_T1_brain.nii', [1 0 0], {'${MNI_REF}1.nii', '${MNI_REF}2.nii', '${MNI_REF}3.nii', '${MNI_REF}4.nii', '${MNI_REF}5.nii', '${MNI_REF}6.nii'}, [5 5 5]); exit;"

      copy_header ${CLASSIFIER}/${SUB}_T1_brain.nii ${CLASSIFIER}/${SUB}_T1.nii

      matlab "-nodesktop -nosplash " <<<"coreg_normalise('${CLASSIFIER}/${SUB}_mean_nds.nii', '${CLASSIFIER}/${SUB}.nii', ${NVOLS}, '${CLASSIFIER}/${SUB}_T1.nii', [0 1 0], {'${MNI_REF}1.nii', '${MNI_REF}2.nii', '${MNI_REF}3.nii', '${MNI_REF}4.nii', '${MNI_REF}5.nii', '${MNI_REF}6.nii'}, [5 5 5]); exit;"

      rm -f ${CLASSIFIER}/c1${SUB}_T1.nii
      rm -f ${CLASSIFIER}/c2${SUB}_T1.nii
      rm -f ${CLASSIFIER}/c3${SUB}_T1.nii
      rm -f ${CLASSIFIER}/c4${SUB}_T1.nii
      rm -f ${CLASSIFIER}/c5${SUB}_T1.nii
      rm -f ${CLASSIFIER}/c6${SUB}_T1.nii



      matlab "-nodesktop -nosplash " <<<"coreg_normalise('${CLASSIFIER}/${SUB}_mean_nds.nii', '${CLASSIFIER}/${SUB}.nii', ${NVOLS}, '${CLASSIFIER}/${SUB}_T1.nii', [0 0 1], {'${MNI_REF}1.nii', '${MNI_REF}2.nii', '${MNI_REF}3.nii', '${MNI_REF}4.nii', '${MNI_REF}5.nii', '${MNI_REF}6.nii'}, [0 0 0]); exit;"

      mv  ${CLASSIFIER}/w${SUB}.nii ${CLASSIFIER}/${SUB}_group.nii
      #3dresample -prefix ${CLASSIFIER}/${SUB}_group.nii -dxyz 2.5 2.5 2.5 -input ${CLASSIFIER}/w${SUB}.nii

      rm -f ${CLASSIFIER}/${SUB}.nii
      rm -f ${CLASSIFIER}/w${SUB}.nii
      rm -f ${CLASSIFIER}/${SUB}_mean_nds.nii
      rm -f ${CLASSIFIER}/${SUB}_mean.nii
      rm -f ${CLASSIFIER}/${SUB}_T1.nii
      rm -f ${CLASSIFIER}/u_rc1${SUB}_T1.nii
      rm -f ${CLASSIFIER}/rc1${SUB}_T1.nii
      rm -f ${CLASSIFIER}/rc2${SUB}_T1.nii
      rm -f ${CLASSIFIER}/rc3${SUB}_T1.nii
      rm -f ${CLASSIFIER}/y_${SUB}_T1.nii
      rm -f ${CLASSIFIER}/iy_${SUB}_T1.nii
      rm -f ${CLASSIFIER}/${SUB}_T1_brain.nii


}
export -f norm_sub

parallel -j10 --line-buffer norm_sub ::: N1_06 N1_14 N1_21 N1_27 N1_34 N1_52 N2_43 N2_05 N2_26 N2_03 N2_40 N2_37 N3_06 N3_15 N3_39 N3_02 N3_18
#B3_17 B2_30 B2_11 B3_37 B3_41 B3_32 B1_07 B3_58 B1_57 B2_48 B1_63 B3_30 B2_60 B1_23 B1_24 B1_09 B2_26
