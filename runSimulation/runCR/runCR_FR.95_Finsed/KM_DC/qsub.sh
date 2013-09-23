Params="0.25 0.3 0.35 0.4 0.483"

for num in $Params;
do
  echo "
  ### Job name
  #PBS -N runCR_KMDC_$num
  ### out files
  #PBS -e ./log/runCR_KMDC_$num.err
  #PBS -o ./log/runCR_KMDC_$num.log
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
  ../../../../winULSAkmeans4b2_DC/winULSAkmeans4b2_DC 195 11 8 0 180 mapFile/mapFile_uni_195_r500/mapFile_uniR500_N195_2.txt $num 1 1 6 0.95 
  " >> ./run_temp.sh
  echo '
  echo End time is `date`
  time2=`date +%s`
  echo Computing time is `echo $time2-$time1 | bc` sec
  ' >> ./run_temp.sh
  qsub run_temp.sh
done

