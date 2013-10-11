
  ### Job name
  #PBS -N ULSA4b2_DC_runHN_FR.95CR.48_N50_25
  ### out files
  #PBS -e ./log/ULSA4b2_DC_runHN_FR.95CR.48_N50_25.err
  #PBS -o ./log/ULSA4b2_DC_runHN_FR.95CR.48_N50_25.log
  ### put the job to which queue (qwork)
  #PBS -q qwork
 
  echo Working directory is $PBS_O_WORKDIR
  cd $PBS_O_WORKDIR
  echo running on host `hostname`
  echo Start time is `date`
  time1=`date +%s`
  echo Directory is `pwd`
 
  ../../../ULSA4b2_DC/ULSA4b2_DC 50 25 8 0 180 mapFile/mapFile_uni_50_r500/mapFile_uniR500_N50_1.txt 0.477 2 1 25 0.95 0 40000 
  

  echo End time is `date`
  time2=`date +%s`
  echo Computing time is `echo $time2-$time1 | bc` sec
  
