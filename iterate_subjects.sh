#!/bin/bash
printf "\033c"




START=$(date -u +%s.%N)

#for i in 69; #`seq $1 $2`;
#do

loop()
{

	#BLOCK=$RANDOM % 9 +1
	
	i=$1
	if (( $i < 10 )); then
		#./proc_rsdata.sh 00$i /home/fsluser/Documents/rs_proc/conf_files/crunch.conf
		SID=00$i
	elif (( $i < 100 )); then
		#./proc_rsdata.sh 0$i /home/fsluser/Documents/rs_proc/conf_files/crunch.conf
		SID=0$i
	else
		#./proc_rsdata.sh $i /home/fsluser/Documents/rs_proc/conf_files/crunch.conf
		SID=$i
	fi
	#./proc_rsdata.sh $i /home/fsluser/Documents/rs_proc/conf_files/crunch.conf

	./proc_rsdata.sh $SID /home/fsluser/Documents/rs_proc/conf_files/crunch_task.conf 1
	
		
	#./proc_rsdata.sh sub$SID /home/fsluser/Documents/rs_proc/conf_files/china.conf
	#./estimate_group_motion.sh sub$SID /home/fsluser/Documents/rs_proc/conf_files/china.conf
	#./estimate_group_motion.sh $SID /home/fsluser/Documents/rs_proc/conf_files/crunch.conf


	#	i=$3
	if [ "$i" -lt "10" ]; then
		ID=0$i
	else
		ID=$i 
	fi

#	if [ "${1}" -eq "B"]; then
#		./estimate_group_motion.sh ${1}${2}_$ID /home/fsluser/Documents/rs_proc/conf_files/rep_impact_belgium.conf
#	else
#		./estimate_group_motion.sh ${1}${2}_$ID /home/fsluser/Documents/rs_proc/conf_files/rep_impact_norway.conf
#	fi

}
export -f loop
#parallel -j3 --line-buffer loop ::: $(seq 1 13) 
parallel -j1 --line-buffer loop :::  2 4 5 6 7 9 10 56 
#parallel -j4 --line-buffer loop ::: B N ::: 1 2 3  ::: $(seq 1 90) 

#done    

END=$(date -u +%s.%N)
DIFF=`echo "(($END - $START))" | bc`
printf "\n\nALL SUBJECTS processed in %.1f seconds\n" $DIFF 
