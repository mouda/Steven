algorithm="MaxEntropy"
epslion="0.001 0.01 0.1 1 10 20"
numNodes="50"
numHeads=(10)
qBits="8"
for token in $algorithm; do 
  j=0
  for eps in $epslion; do
    for (( i = 1; i <= 50; i++ )); do
      ../../../maxEntropy_schedule/simulator -n ${numNodes} -H ${numHeads[$j]} -q $qBits -p 0 -b 180.0 -t 0.1 -I 10.0 -F Kmeans -m \
        ../../mapFile/mapFile_uni_${numNodes}_r150/mapFile_uniR150_N${numNodes}_${i}.txt -f 0.95 -c 0.477 \
        -a $eps -A $token -T 1 -N 1 -E data/GE_${token}_Q${qBits}_N${numNodes}_EPS${eps}_SC.47_TC.${i}.out 
    done
  done
  ((j++));
done
