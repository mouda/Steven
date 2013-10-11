
  ### Job name
  #PBS -N runHN_SADC_23
  ### out files
  #PBS -e ./log/runHN_SADC_23.err
  #PBS -o ./log/runHN_SADC_23.log
  ### put the job to which queue (qwork)
  #PBS -q qwork
 
  echo Working directory is $PBS_O_WORKDIR
  cd $PBS_O_WORKDIR
  echo running on host `hostname`
  echo Start time is `date`
  time1=`date +%s`
  echo Directory is `pwd`
 
  ../../../../ULSA4b2_DC/ULSA4b2_DC 195 23 8 0 180 mapFile/mapFile_uni_195_r500/mapFile_uniR500_N195_2.txt 0.477 2 1 25 0.95 0 40000
  

  echo End time is `date`
  time2=`date +%s`
  echo Computing time is `echo $time2-$time1 | bc` sec
  
