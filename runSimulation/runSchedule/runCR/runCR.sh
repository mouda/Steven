algorithm="Baseline Branchbound GreedyPhysical"
aryCR="0.25 0.3 0.35 0.483"
numNodes="195"
numHeads="39"
qBits="16"
for token in $algorithm; do 
  for CR in $aryCR; do
    for (( i = 1; i <= 50; i++ )); do
      ../../../schedule/simulator -n ${numNodes} -H ${numHeads} -q $qBits -p 0 -b 180.0 -t 0.1 -m \
      ../../mapFile/mapFile_uni_${numNodes}_r150/mapFile_uniR150_N${numNodes}_${i}.txt -f 0.95 -c \
      ${CR} -A $token -T 0.5 -e 0.5 -E data/GE_${token}_Q${qBits}_N${numNodes}_SC.${CR}_TC.5.out -O data/TE_${token}_Q${qBits}_N${numNodes}_SC.${CR}_TC.5.out
    done
  done
done
