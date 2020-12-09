#!/bin/bash
# Clear terminal screen
printf "\033c"

# Change as needed for each project

FIXBIN=/home/luna.kuleuven.be/u0101486/Software/fix/ # PATH where ICA_FIX is installed
WDIR='/home/luna.kuleuven.be/u0101486/workspace/fmri_proc/'
#BASEDIR=/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/tmp/
BASEDIR=/home/luna.kuleuven.be/u0101486/workspace/data/RepImpact/tmp/
#BASEDIR=/home/luna.kuleuven.be/u0101486/workspace/data/RSPET/tmp/
#BASEDIR=/home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/tmp/


#CLASSIFIER=/home/luna.kuleuven.be/u0101486/workspace/data/CRUNCH/fix_classifier_rs/classifier
CLASSIFIER=/home/luna.kuleuven.be/u0101486/workspace/data/RepImpact/fix_classifier_belgium/classifier
#CLASSIFIER=/home/luna.kuleuven.be/u0101486/workspace/data/RepImpact/fix_classifier_norway/classifier
#CLASSIFIER=/home/luna.kuleuven.be/u0101486/workspace/data/RepImpact/fix_classifier_germany_b/classifier
#CLASSIFIER=/home/luna.kuleuven.be/u0101486/workspace/data/RSPET/fix_classifier/classifier
#CLASSIFIER=/home/luna.kuleuven.be/u0101486/workspace/data/CAI_China/fix_classifier_rs/classifier

THR_LIST=(1 3 6 12 20 25 30 35 40 45 50 55 60)
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
#B1_07  OK23
#B1_09  OK23
#B1_18  OK23
#B1_23  EXC [too few components to proper classify]
#B1_24  OK23
#B1_28  OK23
#B1_30  OK23
#B1_55  OK23 [*  Difficult participant]
#B1_63  OK23


#B2_06  OK23
#B2_09  OK23
#B2_11  OK23
#B2_13  OK23
#B2_18  OK23
#B2_24  OK23

#B2_30  OK23 [*]
#B2_53  OK23 [*]
#B2_56  OK23
#B2_61  OK23
#B2_63  EXC [LOTS OF MOV]

#B3_06 OK23
#B3_25 OK23
#B3_27 OK23 [NOT A GOOD SET]
#B3_30 OK23
#B3_32 OK23
#B3_37 OK23
#B3_48 OK23
#B3_53 OK23
#B3_58 OK23
#B3_60 OK23
SUBS=(${BASEDIR}/B1_07/FIX  \
${BASEDIR}/B1_09/FIX  \
${BASEDIR}/B1_18/FIX  \
${BASEDIR}/B1_24/FIX  \
${BASEDIR}/B1_28/FIX  \
${BASEDIR}/B1_30/FIX  \
${BASEDIR}/B1_42/FIX  \
${BASEDIR}/B1_63/FIX  \
${BASEDIR}/B2_06/FIX  \
${BASEDIR}/B2_09/FIX  \
${BASEDIR}/B2_11/FIX  \
${BASEDIR}/B2_13/FIX  \
${BASEDIR}/B2_18/FIX  \
${BASEDIR}/B2_24/FIX  \
${BASEDIR}/B2_56/FIX  \
${BASEDIR}/B3_06/FIX  \
${BASEDIR}/B3_16/FIX  \
${BASEDIR}/B3_30/FIX  \
${BASEDIR}/B3_32/FIX  \
${BASEDIR}/B3_37/FIX  \
${BASEDIR}/B3_53/FIX  \
${BASEDIR}/B3_58/FIX  \
${BASEDIR}/B3_60/FIX  \
${BASEDIR}/B3_62/FIX  )

# RepImpact NORWAY
# N1_02 OK23
# N1_06 OK23
# N1_08 OK23 [High motion]
# N1_12 OK23
# N1_15 OK23
# N1_18 OK23
# N1_21 OK23
# N1_27 OK23
# N1_32 OK23
# N1_36 OK23
# N1_40 OK23
# N1_45 OK23
# N1_52 OK23
# N1_67 OK23
# N1_71 OK23

# N2_03 OK23
# N2_05 OK23
# N2_12 OK23
# N2_21 OK23
# N2_24 OK23
# N2_25 OK23
# N2_37 OK23
# N2_38 OK23 [BAD SUB]
# N2_40 OK23
# N2_42 OK23
# N2_43 OK23
# N2_50 OK23
# N2_62 OK23

# N3_03 OK23
# N3_06 OK23
# N3_09 OK23
# N3_12 OK23 [*]
# N3_14 OK23
# N3_18 OK23
# N3_25 OK23
# N3_26 OK23
# N3_32 OK23
# N3_34 OK23
# N3_40 OK23
# N3_44 OK23 [*]
# N3_58 OK23
# N3_71 OK23
#SUBS=(${BASEDIR}/N1_02/FIX \
#${BASEDIR}/N1_06/FIX \
#{BASEDIR}/N1_08/FIX \
#${BASEDIR}/N1_12/FIX \
#${BASEDIR}/N1_15/FIX \
#${BASEDIR}/N1_18/FIX \
#{BASEDIR}/N1_21/FIX \
#${BASEDIR}/N1_27/FIX \
#{BASEDIR}/N1_32/FIX \
#${BASEDIR}/N1_36/FIX \
#{BASEDIR}/N1_40/FIX \
#${BASEDIR}/N1_45/FIX \
#${BASEDIR}/N1_52/FIX \
#${BASEDIR}/N1_67/FIX \
#${BASEDIR}/N1_71/FIX \
#${BASEDIR}/N2_03/FIX \
#${BASEDIR}/N2_05/FIX \
#${BASEDIR}/N2_12/FIX \
#${BASEDIR}/N2_21/FIX \
#${BASEDIR}/N2_24/FIX \
#${BASEDIR}/N2_25/FIX \
#${BASEDIR}/N2_37/FIX \
#${BASEDIR}/N2_38/FIX \
#${BASEDIR}/N2_40/FIX \
#${BASEDIR}/N2_42/FIX \
#${BASEDIR}/N2_43/FIX \
#${BASEDIR}/N2_50/FIX \
#{BASEDIR}/N2_62/FIX \
#${BASEDIR}/N3_03/FIX \
#${BASEDIR}/N3_06/FIX \
#${BASEDIR}/N3_09/FIX \
#${BASEDIR}/N3_12/FIX \
#${BASEDIR}/N3_14/FIX \
#${BASEDIR}/N3_18/FIX \
#${BASEDIR}/N3_25/FIX \
#${BASEDIR}/N3_26/FIX \
#${BASEDIR}/N3_32/FIX \
#${BASEDIR}/N3_34/FIX \
#{BASEDIR}/N3_40/FIX \
#${BASEDIR}/N3_44/FIX \
#${BASEDIR}/N3_58/FIX \
#${BASEDIR}/N3_71/FIX   )


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
#SUBS=(${BASEDIR}/G1_01/FIX ${BASEDIR}/G1_11/FIX ${BASEDIR}/G1_15/FIX ${BASEDIR}/G1_16/FIX \
#      ${BASEDIR}/G1_30/FIX ${BASEDIR}/G2_01/FIX ${BASEDIR}/G2_12/FIX \
#      ${BASEDIR}/G2_14/FIX ${BASEDIR}/G2_15/FIX ${BASEDIR}/G2_17/FIX ${BASEDIR}/G3_09/FIX \
#      ${BASEDIR}/G3_12/FIX ${BASEDIR}/G3_14/FIX ${BASEDIR}/G3_15/FIX ${BASEDIR}/G3_16/FIX \
#      ${BASEDIR}/G3_17/FIX)

#SUBS=(${BASEDIR}/G1_34/FIX ${BASEDIR}/G1_36/FIX ${BASEDIR}/G1_37/FIX ${BASEDIR}/G1_38/FIX \
      #${BASEDIR}/G1_39/FIX ${BASEDIR}/G1_40/FIX)
# CAI China
#SUBS=(${BASEDIR}/sub101/FIX ${BASEDIR}/sub023/FIX ${BASEDIR}/sub008/FIX ${BASEDIR}/sub115/FIX \
#      ${BASEDIR}/sub122/FIX ${BASEDIR}/sub018/FIX ${BASEDIR}/sub015/FIX ${BASEDIR}/sub110/FIX \
#      ${BASEDIR}/sub109/FIX ${BASEDIR}/sub012/FIX ${BASEDIR}/sub004/FIX ${BASEDIR}/sub124/FIX)


#O01 OK
#O04 Ok
#O06 OK
#O10 OK
#O13 OK
#Y01 OK
#Y02 OK
#Y04 OK
#Y10 OK
#Y14 OK
#Y12 OK
#SUBS=(${BASEDIR}/O01/FIX ${BASEDIR}/O04/FIX ${BASEDIR}/O06/FIX ${BASEDIR}/O10/FIX \
#      ${BASEDIR}/O13/FIX ${BASEDIR}/O15/FIX ${BASEDIR}/Y01/FIX ${BASEDIR}/Y04/FIX ${BASEDIR}/Y10/FIX \
#      ${BASEDIR}/Y12/FIX ${BASEDIR}/Y02/FIX ${BASEDIR}/Y14/FIX ${BASEDIR}/Y07/FIX)
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
