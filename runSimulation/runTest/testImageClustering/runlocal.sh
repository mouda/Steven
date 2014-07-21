tier1Time="90000.0 80000.0 70000.0 60000.0 50000.0"
for token in $tier1Time
do
  ../../../cluster/genCluster -n 24 -H 4 -q 32 -p 30 -b 180.0 -t 20.0 -I ${token} -m  \
    paper720_30cam_pos.txt -f 1 -c 0.477  -T  1.0 -F  ImageSource -N 5  -C CS_${token}.out
done
