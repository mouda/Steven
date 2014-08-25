algorithm="Branchbound GreedyPhysical SumRate"
numNodes="50 100 150 195"
numHeads=(10 20 30 39)
qBits="8"
for token in $algorithm; do 
  j=0
  for num in $numNodes; do
    for (( i = 1; i <= 50; i++ )); do
      ../../../maxEntropy_schedule/simulator -n ${num} -H ${numHeads[$j]} -q $qBits -p 0 -b 180.0 -t 0.1 -I 10.0 -F Kmeans -m \
      ../../mapFile/mapFile_uni_${num}_r150/mapFile_uniR150_N${num}_${i}.txt -f 0.95 -c 0.477 \
      -a 0.001 -A $token -T 0.5 -N 1 -E data/GE_${token}_Q${qBits}_N${num}_SC.47_TC.1_Idx${i}.out 
    done
    ((j++));
  done
done
