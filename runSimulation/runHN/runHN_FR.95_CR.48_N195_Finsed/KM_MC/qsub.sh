for ((addue=3;addue<=23;addue=addue+2))
do
  echo "
  ### Job name
  #PBS -N runHN_KMMC_$addue
  ### out files
  #PBS -e ./log/runHN_KMMC_$addue.err
  #PBS -o ./log/runHN_KMMC_$addue.log
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
  ../../../../winULSAkmeans2i_MC/winULSAkmeans2i_MC 195 $addue 8 0 180 mapFile/mapFile_uni_195_r500/mapFile_uniR500_N195_2.txt 0.477 2 1 6 0.95 
  " >> ./run_temp.sh
  echo '
  echo End time is `date`
  time2=`date +%s`
  echo Computing time is `echo $time2-$time1 | bc` sec
  ' >> ./run_temp.sh
  qsub run_temp.sh
done

