algorithm="Baseline Branchbound GreedyPhysical"
qBits="16"
for token in $algorithm; do
  for (( i = 1; i <= 50; i++ )); do
    ../../../schedule/simulator -n 50 -H 10 -q $qBits -p 0 -b 180.0 -t 0.1 -m \
    ../../mapFile/mapFile_uni_50_r150/mapFile_uniR150_N50_${i}.txt -f 0.95 -c \
    0.477 -A $token -T 0.5 -e 0.5 -E entropy_${token}_Q_${qBits}_N50.out -M MSE_${token}_N50.out
  done
done
