algorithm="MaxSNR"
numNodes="50 100 150 195"
numHeads=(10 20 30 39)
qBits="8"
for token in $algorithm; do 
  j=0
  for num in $numNodes; do
    for (( i = 1; i <= 50; i++ )); do
      echo "=========== ${token}, Nudes: ${num}, Index : $i ========="
      echo "
      ### Job name
      #PBS -N simulator_${token}_Q${qBits}_N${num}_SC.47_TC.1_Idx${i} 
      ### out files
      #PBS -e ./log/simulator_${token}_Q${qBits}_N${num}_SC.47_TC.1_Idx${i}.err
      #PBS -o ./log/simulator_${token}_Q${qBits}_N${num}_SC.47_TC.1_Idx${i}.log
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
      ../../../maxEntropy_schedule/simulator -n ${num} -H ${numHeads[$j]} -q $qBits -p 0 -b 180.0 -t 0.1 -I 10.0 -F Kmeans -m \
      ../../mapFile/mapFile_uni_${num}_r150/mapFile_uniR150_N${num}_${i}.txt -f 0.95 -c 0.477 \
      -a 0.001 -A $token -T 0.5 -N 1 -E data/GE_${token}_Q${qBits}_N${num}_SC.47_TC.1_Idx${i}.out 
      " >> ./run_temp.sh
      echo '
      echo End time is `date`
      time2=`date +%s`
      echo Computing time is `echo $time2-$time1 | bc` sec
      ' >> ./run_temp.sh
      qsub run_temp.sh
    done
    ((j++));
  done
done
