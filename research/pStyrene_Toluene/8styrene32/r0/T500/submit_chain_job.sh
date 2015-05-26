#!/bin/sh

istart=2 #Start run number istart
iend=10 #End run number
for i in `seq $istart $iend`;do
 /bin/cp equil.I.in equil.${i}.in
 perl -p -i -e "s/_I_/$i/g" equil.${i}.in
 /bin/cp equil_hopper_I.pbs equil_hopper_${i}.pbs
 perl -p -i -e "s/_I_/$i/g" equil_hopper_${i}.pbs
 if [ $i = $istart ];then
  PBS_JOBID=`qsub -q regular equil_hopper_${i}.pbs`
  #echo 'no'
 else
  PBS_JOBID=`qsub -q regular -W depend=afterok:$PBS_JOBID equil_hopper_${i}.pbs`
  #echo 'yes'
 fi
 echo $PBS_JOBID
 sleep 1s
done
