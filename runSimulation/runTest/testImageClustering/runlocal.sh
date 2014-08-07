tier1Time="90000.0 80000.0 70000.0 60000.0 50000.0"
for token in $tier1Time
do
  ../../../cluster/genCluster -n 24 -H 4 -q 32 -p 30 -b 180.0 -t 20.0 -I ${token} -m  \
    paper720_30cam_pos.txt -f 1 -c 0.477  -T  1.0 -F  ImageSource -N 5  -C CS_${token}.out
done
../../../image_cluster/genCluster -n 100 -H 17  -p 30 -b 180 -t 170000.0 -I 21000.0 -m source/100cam_r2000_map.out -f 1.0 -c source/100cam_r2000_corr.out  -i source/100cam_r2000_idt.out -F  ImageBaseline -N 5  -C CS_N100_R2000_1.0_H60_N5_baseline.out
../../../image_cluster/genCluster -n 100 -H 60  -p 30 -b 180 -t 170000.0 -I 21000.0 -m source/100cam_r2000_map.out -f 1.0 -c source/100cam_r2000_corr.out  -i source/100cam_r2000_idt.out -F  ImageSource -N 5  -C CS_test_CSA.out
