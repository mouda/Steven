algorithm="Branchbound"
numNodes="195"
numHeads=(39)
qBits="16"
for token in $algorithm; do 
  j=0
  for num in $numNodes; do
    echo ${numHeads[$j]}
    for (( i = 1; i <= 50; i++ )); do
      ../../../schedule/simulator -n ${num} -H ${numHeads[$j]} -q $qBits -p 0 -b 180.0 -t 0.1 -m \
      ../../mapFile/mapFile_uni_${num}_r150/mapFile_uniR150_N${num}_${i}.txt -f 0.95 -c \
      0.477 -A $token -T 0.5 -e 0.5 -E data/entropy_${token}_Q_${qBits}_N${num}.out -M data/MSE_${token}_N${num}.out
    done
    ((j++));
  done
done
