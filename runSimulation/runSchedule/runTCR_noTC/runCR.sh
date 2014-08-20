algorithm="Baseline Branchbound GreedyPhysical"
SCR="0.483"
aryTCR="0.001"
#aryTCR="0.01 0.03 0.1 0.32 1 3.16 10"
#aryTCR="0.01 0.02 0.03 0.06 0.1 0.32 1 3.16 10"
numNodes="195"
numHeads="39"
qBits="16"
for token in $algorithm; do 
  for TCR in $aryTCR; do
    for (( i = 1; i <= 50; i++ )); do
      ../../../schedule/simulator -n ${numNodes} -H ${numHeads} -q $qBits -p 0 -b 180.0 -t 0.1 -m \
      ../../mapFile/mapFile_uni_${numNodes}_r150/mapFile_uniR150_N${numNodes}_${i}.txt -f 0.95 -c \
      ${SCR} -A $token -T $TCR -e 0.5 -E data/GE_${token}_Q${qBits}_N${numNodes}_SC.${SCR}_TC.${TCR}.out -O data/TE_${token}_Q${qBits}_N${numNodes}_SC.${SCR}_TC.${TCR}.out
    done
  done
done
