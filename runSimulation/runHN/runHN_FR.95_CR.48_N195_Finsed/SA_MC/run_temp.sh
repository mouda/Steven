
  ### Job name
  #PBS -N winULSAkmeans2i_MC_runHN_FR.95CR.48RanMore__Finsed_25
  ### out files
  #PBS -e ./log/winULSA2i2_MC_runHN_FR.95CR.48RanMore__Finsed_25.err
  #PBS -o ./log/winULSA2i2_MC_runHN_FR.95CR.48RanMore__Finsed_25.log
  ### put the job to which queue (qwork)
  #PBS -q qwork
 
  echo Working directory is $PBS_O_WORKDIR
  cd $PBS_O_WORKDIR
  echo running on host `hostname`
  echo Start time is `date`
  time1=`date +%s`
  echo Directory is `pwd`
 
  ./winULSA2i2_MC 195 25 8 0 180 mapFile/mapFile_uni_195_r500/mapFile_uniR500_N195_2.txt 0.477 2 1 6 0.95
  

  echo End time is `date`
  time2=`date +%s`
  echo Computing time is `echo $time2-$time1 | bc` sec
  
