#!/bin/bash
# Clear terminal screen
printf "\033c"

# Change as needed for each project

FIXBIN=/home/luna.kuleuven.be/u0101486/Software/fix/ # PATH where ICA_FIX is installed
WDIR='/home/luna.kuleuven.be/u0101486/workspace/fmri_proc/'
#BASEDIR=/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/
BASEDIR=/home/luna.kuleuven.be/u0101486/workspace/data/RepImpact/tmp/
#BASEDIR=/home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/tmp/


#CLASSIFIER=/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/fix_classifier_rs/classifier
#CLASSIFIER=/home/luna.kuleuven.be/u0101486/workspace/data/RepImpact/fix_classifier_belgium_b/classifier
#CLASSIFIER=/home/luna.kuleuven.be/u0101486/workspace/data/RepImpact/fix_classifier_norway_b/classifier
CLASSIFIER=/home/luna.kuleuven.be/u0101486/workspace/data/RepImpact/fix_classifier_germany_b/classifier
#CLASSIFIER=/home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/fix_classifier_rs/classifier

THR_LIST=(1 5 10 15 25 35 50 60)
#1 - RS005 OK
#2 - RS006 OK
#3 - RS010 OK
#4 - RS015 OK
#5 - RS016 OK
#6 - RS018 OK
#7 - RS020 OK
#8 - RS023 Ok
#9 - RS026 OK
#10 - RS027 OK
#11 - RS029 OK
#12 - RS032 OK
#13 - RS036 OK
#14 - RS050 OK
#15 - RS061 OK
#16 - RS075 OK
#17 - RS081 OK
#18 - RS086 OK
#19 - RS087 OK
#20 - RS103
# CRUNCH
#SUBS=(${BASEDIR}/RS005/FIX ${BASEDIR}/RS006/FIX ${BASEDIR}/RS010/FIX ${BASEDIR}/RS015/FIX \
#      ${BASEDIR}/RS016/FIX ${BASEDIR}/RS018/FIX ${BASEDIR}/RS020/FIX ${BASEDIR}/RS023/FIX \
#      ${BASEDIR}/RS026/FIX ${BASEDIR}/RS027/FIX ${BASEDIR}/RS029/FIX ${BASEDIR}/RS032/FIX \
#      ${BASEDIR}/RS036/FIX ${BASEDIR}/RS050/FIX ${BASEDIR}/RS061/FIX ${BASEDIR}/RS075/FIX \
#      ${BASEDIR}/RS081/FIX ${BASEDIR}/RS086/FIX ${BASEDIR}/RS087/FIX ${BASEDIR}/RS103/FIX)



# RepImpact BELGIUM
#B1_07  OK
#B1_09  OK
#B1_23  OK
#B1_24  OK
#B1_28  OK
#B1_30  OK
#B1_63  OK

#B2_09  OK
#B2_11  OK

#B2_24  OK
#B2_26  OK
#B2_30  OK
#B2_50   --> Does not exist?
#B2_53  OK
#B2_56  OK

#B3_06 OK
#B3_14 OK
#B3_30 OK
#B3_32 OK
#B3_37 OK
#B3_41 OK
#B3_58 OK
#SUBS=(${BASEDIR}/B1_07/FIX ${BASEDIR}/B1_09/FIX ${BASEDIR}/B1_23/FIX ${BASEDIR}/B1_24/FIX \
#      ${BASEDIR}/B1_28/FIX ${BASEDIR}/B1_30/FIX ${BASEDIR}/B1_55/FIX ${BASEDIR}/B1_63/FIX \
#      ${BASEDIR}/B2_09/FIX ${BASEDIR}/B2_11/FIX ${BASEDIR}/B2_13/FIX ${BASEDIR}/B2_24/FIX   \
#      ${BASEDIR}/B2_30/FIX ${BASEDIR}/B2_53/FIX ${BASEDIR}/B2_56/FIX ${BASEDIR}/B2_63/FIX \
#      ${BASEDIR}/B3_06/FIX ${BASEDIR}/B3_27/FIX ${BASEDIR}/B3_30/FIX ${BASEDIR}/B3_32/FIX   \
#      ${BASEDIR}/B3_37/FIX ${BASEDIR}/B3_53/FIX ${BASEDIR}/B3_58/FIX ${BASEDIR}/B3_60/FIX )

# RepImpact NORWAY
# N1_02 OK
# N1_06 OK
# N1_08 OK
# N1_21 OK
# N1_27 OK
# N1_36 OK
# N1_52 OK
# N2_03 OK
# N2_05 OK
# N2_12 OK
# N2_25 OK

# N2_37 OK
# N2_40 OK
# N2_42 OK
# N2_43 OK
# N2_62 OK

# N3_06 OK
# N3_12 OK
# N3_25 OK
# N3_32 OK

# N3_34 OK
# N3_40 OK
# N3_58 OK
#SUBS=(${BASEDIR}/N1_02/FIX ${BASEDIR}/N1_06/FIX ${BASEDIR}/N1_08/FIX ${BASEDIR}/N1_12/FIX \
#      ${BASEDIR}/N1_18/FIX  ${BASEDIR}/N1_21/FIX \
#      ${BASEDIR}/N1_27/FIX ${BASEDIR}/N1_36/FIX ${BASEDIR}/N1_52/FIX \
#      ${BASEDIR}/N1_67/FIX ${BASEDIR}/N2_03/FIX \
#      ${BASEDIR}/N2_05/FIX ${BASEDIR}/N2_12/FIX ${BASEDIR}/N2_21/FIX \
#      ${BASEDIR}/N2_25/FIX ${BASEDIR}/N2_38/FIX ${BASEDIR}/N2_37/FIX ${BASEDIR}/N2_40/FIX \
#      ${BASEDIR}/N2_42/FIX ${BASEDIR}/N2_43/FIX ${BASEDIR}/N2_62/FIX \
#      ${BASEDIR}/N3_06/FIX ${BASEDIR}/N3_09/FIX ${BASEDIR}/N3_12/FIX \
#      ${BASEDIR}/N3_25/FIX ${BASEDIR}/N3_26/FIX ${BASEDIR}/N3_32/FIX \
#      ${BASEDIR}/N3_34/FIX ${BASEDIR}/N3_40/FIX  ${BASEDIR}/N3_44/FIX ${BASEDIR}/N3_58/FIX \
#      ${BASEDIR}/N3_71/FIX )


# G1_01 OK2
# G1_11 OK2
# G1_15 OK2
# G1_16 OK2
# G1_23 OK2
# G1_30 OK2


# G2_01 OK2
# G2_08 OK2
# G2_12 OK2
# G2_14 OK2
# G2_15 OK2
# G2_17 OK2

# G3_09 OK2
# G3_12 OK2
# G3_14 OK2
# G3_15 OK2
# G3_16 OK2
# G3_17 OK2
SUBS=(${BASEDIR}/G1_01/FIX ${BASEDIR}/G1_11/FIX ${BASEDIR}/G1_15/FIX ${BASEDIR}/G1_16/FIX \
      ${BASEDIR}/G1_30/FIX ${BASEDIR}/G2_01/FIX ${BASEDIR}/G2_12/FIX \
      ${BASEDIR}/G2_14/FIX ${BASEDIR}/G2_15/FIX ${BASEDIR}/G2_17/FIX ${BASEDIR}/G3_09/FIX \
      ${BASEDIR}/G3_12/FIX ${BASEDIR}/G3_14/FIX ${BASEDIR}/G3_15/FIX ${BASEDIR}/G3_16/FIX \
      ${BASEDIR}/G3_17/FIX)
# CAI China
#SUBS=(${BASEDIR}/sub101/FIX ${BASEDIR}/sub023/FIX ${BASEDIR}/sub008/FIX ${BASEDIR}/sub115/FIX \
#      ${BASEDIR}/sub122/FIX ${BASEDIR}/sub018/FIX ${BASEDIR}/sub015/FIX ${BASEDIR}/sub110/FIX \
#      ${BASEDIR}/sub109/FIX ${BASEDIR}/sub012/FIX ${BASEDIR}/sub004/FIX ${BASEDIR}/sub124/FIX)
NSUBS=0


for SUB in ${SUBS[@]}; do
  NSUBS=$((NSUBS+1))
done

# LEAVE ONE OUT PROCEDURE
for ONE_LEFT in $(seq 1 $NSUBS); do

      ONE_LEFT=$((ONE_LEFT-1))

      rm -f ${CLASSIFIER}.RData
      rm -r -f ${CLASSIFIER}

      FIX_CMD="${FIXBIN}/fix -t ${CLASSIFIER}"

      SI=0

      for SUB in ${SUBS[@]}; do
            if [ "${ONE_LEFT}" -eq "${SI}" ]; then
                  echo "Leaving ${SUB} out"
                  OUT_SUB=${SUB}
            else
                  FIX_CMD=`echo "${FIX_CMD} " ${SUB} " "`
            fi
            SI=$((SI+1))
      done

      eval ${FIX_CMD}

      for THR in ${THR_LIST[@]}; do
        eval "${FIXBIN}/fix -c ${OUT_SUB} ${CLASSIFIER}.RData ${THR}"
        python ${WDIR}/compare_fix_class.py -hand ${OUT_SUB}/hand_labels_noise.txt -fix ${OUT_SUB}/fix4melview_classifier_thr${THR}.txt \
                -out ${CLASSIFIER}_${ONE_LEFT}_${THR}.txt -thr ${THR}
        rm ${OUT_SUB}/fix4melview_classifier_thr${THR}.txt
      done
done


#To calculate overall LOO performance, run the loo_performance.py script


# Train the classifier with all subjects
rm -f ${CLASSIFIER}.RData
rm -r -f ${CLASSIFIER}
FIX_CMD="${FIXBIN}/fix -t ${CLASSIFIER}"
for SUB in ${SUBS[@]}; do
      FIX_CMD=`echo "${FIX_CMD} " ${SUB} " "`
done
eval ${FIX_CMD}
