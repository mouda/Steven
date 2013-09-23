Params="0.55 0.65 0.75 0.85 0.95" # Fidility Ration

for num in $Params;
do
  echo "
  ### Job name
  #PBS -N runFR_KMMC_$num
  ### out files
  #PBS -e ./log/runFR_KMMC_$num.err
  #PBS -o ./log/runFR_KMMC_$num.log
  ### put the job to which queue (qwork)
  #PBS -q qwork" > ./run_temp.sh
  echo ' 
  echo Working directory is $PBS_O_WORKDIR
  cd $PBS_O_WORKDIR
  echo running on host `hostname`
  echo Start time is `date`
  time1=`date +%s`
  echo Directory is `pwd`' >> ./run_temp.sh
  echo " 
  ../../../../winULSAkmeans2i_MC/winULSAkmeans2i_MC 195 5 8 0 180 mapFile/mapFile_uni_195_r500/mapFile_uniR500_N195_2.txt 0.48 2 1 6 $num 
  " >> ./run_temp.sh
  echo '
  echo End time is `date`
  time2=`date +%s`
  echo Computing time is `echo $time2-$time1 | bc` sec
  ' >> ./run_temp.sh
  qsub run_temp.sh
done

