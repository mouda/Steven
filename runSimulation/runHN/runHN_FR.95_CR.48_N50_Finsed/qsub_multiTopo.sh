rm tmpAll/*

for ((addue=3;addue<=23;addue=addue+2))
do
  for ((topoIdx=1;topoIdx<=50;topoIdx=topoIdx+1))
  do
    echo "
    ### Job name
    #PBS -N ULSA4b2_DC_runHN_N50_$addue
    ### out files
    #PBS -e ./log/ULSA4b2_DC_runHN_N50_${addue}_${topoIdx}.err
    #PBS -o ./log/ULSA4b2_DC_runHN_N50_${addue}_${topoIdx}.log
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
    ../../../ULSA4b2_DC/ULSA4b2_DC 50 $addue 8 0 180 ../../mapFile/mapFile_uni_50_r500/mapFile_uniR500_N50_${topoIdx}.txt 0.477 2 1 1 0.95 0 20000 
    " >> ./run_temp.sh
    echo '
    echo End time is `date`
    time2=`date +%s`
    echo Computing time is `echo $time2-$time1 | bc` sec
    ' >> ./run_temp.sh
    qsub run_temp.sh
  done
done

# Proccess remain files

cat tmpAll/ULSA4b2_All_N50_BW180.0PW0.001_FR0.95_r500.0.* > ULSA4b2_All_N50_BW180.0PW0.001_FR0.95_r500.0.txt 


