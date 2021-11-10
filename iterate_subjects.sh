#!/bin/bash
printf "\033c"
START=$(date -u +%s.%N)


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
	./template_generation/import_t1s.sh ${1}${ID}  /media/thiago/EXTRALINUX/fmri_proc/import_scripts/connect_ex_t1.sh /media/thiago/EXTRALINUX/data/ConnectEx/Template/anat


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
  ./proc_main.sh ${1}${ID}  /media/thiago/EXTRALINUX/fmri_proc/conf_files/connectex.conf

  #./proc_main.sh ${1}${ID}_noPhys  /home/luna.kuleuven.be/u0101486/workspace/fmri_proc/conf_files/connectex_nophys.conf


}
export -f loop_connex


rm -f /home/thiago/.afni.log

parallel -j3 --line-buffer loop_connex ::: A ::: 90


#parallel -j10 --line-buffer loop_t1_connex ::: A Y ::: $(seq 1 120)
#parallel -j2 --line-buffer loop_t1 ::: N  ::: 1 2 3 ::: $(seq 1 75)
#matlab "-nodesktop -nosplash " <<<"generate_dartel_template('/home/luna.kuleuven.be/u0101486/workspace/data/RepImpact/Template/anat/'); exit;"

END=$(date -u +%s.%N)
DIFF=`echo "(($END - $START))" | bc | awk '{printf "All subs %f s", $0}'`

#printf "\n\nALL SUBJECTS processed in %.1f seconds\n" $DIFF
