#!/bin/bash
printf "\033c"




START=$(date -u +%s.%N)

loop_pet()
{
	i=$2
	if [ "$i" -lt "10" ]; then
		ID=0$i
	else
		ID=$i
	fi

	./proc_main.sh ${1}${ID}  /home/luna.kuleuven.be/u0101486/workspace/fmri_proc/conf_files/rs_pet.conf
}
export -f loop_pet


loop_t1()
{
	i=$3
	if [ "$i" -lt "10" ]; then
		ID=0$i
	else
		ID=$i
	fi

  echo "${1}${2}_$ID"
	if test "${1}" == "B"; then
		./import_t1s.sh ${1}${2}_${ID}  /home/luna.kuleuven.be/u0101486/workspace/fmri_proc/import_scripts/repimpact_t1.sh /home/luna.kuleuven.be/u0101486/workspace/data/RepImpact/Template/anat/
	fi

  if test "${1}" == "N"; then
		./import_t1s.sh ${1}${2}_${ID}  /home/luna.kuleuven.be/u0101486/workspace/fmri_proc/import_scripts/repimpact_norway_t1.sh /home/luna.kuleuven.be/u0101486/workspace/data/RepImpact/Template/anat/
	fi

  if test "${1}" == "G"; then
		./import_t1s.sh ${1}${2}_${ID}  /home/luna.kuleuven.be/u0101486/workspace/fmri_proc/import_scripts/repimpact_germany_t1.sh /home/luna.kuleuven.be/u0101486/workspace/data/RepImpact/Template/anat/
	fi

}
export -f loop_t1

loop_t1_connex()
{
	i=$2

  if [ "$i" -lt "10" ]; then
    ID=00$i
  elif [ "$i" -lt "100" ]; then
		ID=0$i
	else
		ID=$i
	fi

  echo "${1}${2}_$ID"
	./template_generation/import_t1s.sh ${1}${ID}  /home/luna.kuleuven.be/u0101486/workspace/fmri_proc/import_scripts/connect_ex_t1.sh /home/luna.kuleuven.be/u0101486/workspace/data/ConnectEx/Template/anat


}
export -f loop_t1_connex


loop_connex()
{
	i=$2
  P=$3
  if [ "$i" -lt "10" ]; then
    ID=00$i
  elif [ "$i" -lt "100" ]; then
		ID=0$i
	else
		ID=$i
	fi



  echo "${1}$ID"
  ./proc_main.sh ${1}${ID}  /home/luna.kuleuven.be/u0101486/workspace/fmri_proc/conf_files/connectex.conf

  #./proc_main.sh ${1}${ID}_noPhys  /home/luna.kuleuven.be/u0101486/workspace/fmri_proc/conf_files/connectex_nophys.conf


}
export -f loop_connex

loop()
{
	i=$3
	if [ "$i" -lt "10" ]; then
		ID=0$i
	else
		ID=$i
	fi

  #echo "${1}${2}_$ID"
	if test "${1}" == "B"; then
		./proc_main.sh ${1}${2}_${ID}  /home/luna.kuleuven.be/u0101486/workspace/fmri_proc/conf_files/rep_impact_belgium.conf
	fi

  if test "${1}" == "N"; then
		./proc_main.sh ${1}${2}_${ID}  /home/luna.kuleuven.be/u0101486/workspace/fmri_proc/conf_files/rep_impact_norway.conf
	fi

  if test "${1}" == "G"; then
    if [ "$i" -lt "32" ]; then
		    ./proc_main.sh ${1}${2}_${ID}  /home/luna.kuleuven.be/u0101486/workspace/fmri_proc/conf_files/rep_impact_germany.conf
    else
        ./proc_main.sh ${1}${2}_${ID}  /home/luna.kuleuven.be/u0101486/workspace/fmri_proc/conf_files/rep_impact_germany_21.conf
    fi
	fi

}
export -f loop

loop_crunch()
{

	#BLOCK=$RANDOM % 9 +1

	i=$1
	if (( $i < 10 )); then
		SID=00$i
	elif (( $i < 100 )); then
		SID=0$i
	else
		SID=$i
	fi
	./proc_main.sh $SID /home/luna.kuleuven.be/u0101486/workspace/fmri_proc/conf_files/crunch.conf 3

}
export -f loop_crunch

loop_cai_china()
{

	#BLOCK=$RANDOM % 9 +1
  GROUP=$1
  i=$2

  if [ "${GROUP}" -eq "0" ]; then
      if [ "${i}" -eq "4" ] || [ "${i}" -eq "10" ] || [ "${i}" -eq "17" ] || \
      [ "${i}" -eq "21" ] || [ "${i}" -eq "22" ]; then
        exit
      fi
  fi

  if [ "${GROUP}" -eq "1" ]; then
      if [ "${i}" -eq "4" ] || [ "${i}" -eq "10" ] || [ "${i}" -eq "13" ]; then
        exit
      fi
  fi

	if (( $i < 10 )); then
		SID=sub${GROUP}0$i
	else
		SID=sub${GROUP}$i
	fi
	./proc_main.sh $SID /home/luna.kuleuven.be/u0101486/workspace/fmri_proc/conf_files/china.conf 1

}
export -f loop_cai_china

#$(seq 1 106)
#parallel -j10 --line-buffer loop_crunch ::: $(seq 1 106) #$(seq 1 106)

#matlab "-nodesktop -nosplash " <<<"run_icap; exit;"

rm -f /home/luna.kuleuven.be/u0101486/.afni.log

#parallel -j10 --line-buffer loop ::: G B N ::: 1 2 3 ::: $(seq 1 75)
#parallel -j10 --line-buffer loop ::: B N ::: 1 2 3 ::: $(seq 1 75)
#parallel -j8 --line-buffer loop_cai_china ::: 0 ::: 1
#parallel -j8 --line-buffer loop_cai_china ::: 0 ::: 16
#parallel -j3 --line-buffer loop_crunch ::: 20

#parallel -j3 --line-buffer loop_crunch ::: 5 6 7 8 9 10 11 12 13  #$(seq 1 106)

# O Y
parallel -j5 --line-buffer loop_pet ::: O Y ::: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18
#parallel -j5 --line-buffer loop_pet ::: O ::: 2


#parallel -j8 --line-buffer loop ::: B  ::: 3 :::  24
#parallel -j3 --line-buffer loop_crunch ::: $(seq 1 106)

#N1_06 N1_14 N1_21 N1_27 \
#      N1_34 N1_52 N2_43 N2_05 \
#      N2_26 N2_03 N2_40 N2_37 \
#      N3_06 N3_15 N3_39 N3_02 \
#      N3_18


#parallel -j12 --line-buffer loop ::: B N G ::: 1 2 3 ::: $(seq 1 75)
#parallel -j12 --line-buffer loop ::: G ::: 2 ::: 15

#parallel -j10 --line-buffer loop ::: B N G ::: 1 2 3 ::: $(seq 1 75)
#parallel -j12 --line-buffer loop ::: B ::: 1 ::: 7
#parallel -j12 --line-buffer loop ::: B  ::: 2 3 ::: $(seq 1 75)
##parallel -j5 --line-buffer loop ::: B  :::  3 ::: $(seq 20 75)
#parallel -j12 --line-buffer loop ::: N ::: 1  ::: 6 21 27 34 52
#parallel -j12 --line-buffer loop ::: N ::: 2  ::: 43 5 26 3 40 37
#parallel -j12 --line-buffer loop ::: N ::: 3  ::: 6 15 39 2 18
#parallel -j12 --line-buffer loop ::: B ::: 1 2 3 ::: $(seq 1 75)
#parallel -j8 --line-buffer loop ::: B ::: 1 ::: 2
#parallel -j8 --line-buffer loop ::: B ::: 1 ::: 23 24 28 57 63
#parallel -j8 --line-buffer loop ::: B ::: 2 ::: 9 11 26 30 48 60
#parallel -j8 --line-buffer loop ::: B ::: 3 ::: 17 37 41 30 32 58

#parallel -j3 --line-buffer loop_connex ::: A Y ::: $(seq 1 3)
# PHYSIO is failing... see why
#parallel -j3 --line-buffer loop_connex ::: A ::: $(seq 1 50)


#parallel -j10 --line-buffer loop_t1_connex ::: A Y ::: $(seq 1 120)
#parallel -j8 --line-buffer loop_t1 ::: G  ::: 1 2 3 ::: $(seq 1 75)
#matlab "-nodesktop -nosplash " <<<"generate_dartel_template('/home/luna.kuleuven.be/u0101486/workspace/data/RepImpact/Template/anat/'); exit;"


#parallel -j10 --line-buffer loop_crunch ::: $(seq 1 106) #14 16 17 36 41 46 51 52 60 62 63 70 85 97 102
#parallel -j3 --line-buffer loop ::: N ::: 1 2 3 ::: $(seq 1 75)
#parallel -j2 --line-buffer loop_cai_china ::: 1  ::: 1
#parallel -j6 --line-buffer loop_cai_china ::: 0 1 ::: $(seq 1 30)

#parallel -j10 --line-buffer loop_cai_china ::: 0 1 ::: $(seq 1 30)
#parallel -j10 --line-buffer loop_cai_china ::: 0 ::: 6 8 9
#parallel -j10 --line-buffer loop_cai_china ::: 1 ::: 3 4 10 15
#parallel -j10 --line-buffer loop_cai_china ::: 0 ::: 8
#parallel -j10 --line-buffer loop_cai_china ::: 0 ::: 9
#parallel -j10 --line-buffer loop_cai_china ::: 0 ::: 1





#done

END=$(date -u +%s.%N)
DIFF=`echo "(($END - $START))" | bc`
printf "\n\nALL SUBJECTS processed in %.1f seconds\n" $DIFF
