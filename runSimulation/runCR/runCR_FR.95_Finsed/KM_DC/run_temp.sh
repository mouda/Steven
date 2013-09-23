
  ### Job name
  #PBS -N runCR_KMMC_0.483
  ### out files
  #PBS -e ./log/runCR_KMMC_0.483.err
  #PBS -o ./log/runCR_KMMC_0.483.log
  ### put the job to which queue (qwork)
  #PBS -q qwork
 
  echo Working directory is $PBS_O_WORKDIR
  cd $PBS_O_WORKDIR
  echo running on host `hostname`
  echo Start time is `date`
  time1=`date +%s`
  echo Directory is `pwd`
 
  ../../../../winULSAkmeans4b2_DC/winULSAkmeans4b2_DC 195 11 8 0 180 mapFile/mapFile_uni_195_r500/mapFile_uniR500_N195_2.txt 0.483 1 1 6 0.95 
  

  echo End time is `date`
  time2=`date +%s`
  echo Computing time is `echo $time2-$time1 | bc` sec
  