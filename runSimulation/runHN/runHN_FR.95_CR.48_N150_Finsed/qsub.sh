for ((addue=3;addue<=23;addue=addue+2))
do
  echo "
  ### Job name
  #PBS -N runHN_SADC_N150_$addue
  ### out files
  #PBS -e ./log/ULSA4b2_DC_runHN_FR.95CR.48_N150_$addue.err
  #PBS -o ./log/ULSA4b2_DC_runHN_FR.95CR.48_N150_$addue.log
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
  ../../../ULSA4b2_DC/ULSA4b2_DC 150 $addue 8 0 180 mapFile/mapFile_uni_150_r500/mapFile_uniR500_N150_1.txt 0.477 2 1 12 0.95 0 80000 
  " >> ./run_temp.sh
  echo '
  echo End time is `date`
  time2=`date +%s`
  echo Computing time is `echo $time2-$time1 | bc` sec
  ' >> ./run_temp.sh
  qsub run_temp.sh
done

