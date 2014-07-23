tier1Time="90000.0 80000.0 70000.0 60000.0 50000.0"
totalTime="100000.0"
tier2slot="5.0"

for token in $tier1Time
do
  for (( i = 1; i <= 100; i++ )); do
    tier2slotTime=`echo "($totalTime-$token)/$tier2slot" | bc`
    ../../../cluster/genCluster -n 24 -H 4 -q 32 -p 30 -b 180.0 -t ${tier2slotTime} -I ${token} -m  \
      paper720_30cam_pos.txt -f 1 -c 0.477  -T  1.0 -F  ImageBaseline -N 5  -C data/CS_baseline_${token}_${i}.out
  done
done
