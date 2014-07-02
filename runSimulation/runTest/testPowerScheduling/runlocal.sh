#!/bin/bash

#for (( i = 1; i < 6; i++ )); do
#  echo "========== # of slots : $i ==========="
#  echo "========== # of slots : $i ===========" >> time_record.out
#  (time ../../../schedule/simulator -n 50 -H 5 -q 32 -p 0 -b 180.0 -t 0.1 -m   ../../mapFile/mapFile_uni_50_r150/mapFile_uniR150_N50_1.txt -f 0.95 -c 0.477 -A MinPower -T  1.0 -F Kmeans -N $i  -S solution_N50_Slot$i.out -P power_N50_Slot$i.out) 2>> time_record.out
#done

#for (( i = 1; i < 6; i++ )); do
#  echo "========== # of slots : $i ==========="
#  echo "========== # of slots : $i ===========" >> time_record.out
#  (time ../../../schedule/simulator -n 195 -H 10 -q 32 -p 0 -b 180.0 -t 0.1 -m   ../../mapFile/mapFile_uni_195_r150/mapFile_uniR150_N195_1.txt -f 0.95 -c 0.477 -A MinPower -T  1.0 -F Kmeans -N $i  -S solution_N195.out -P power_N195.out) 2>> time_record_N195.out
#done

#for (( i = 1; i < 6; i++ )); do
#  echo "========== # of slots : $i ==========="
#  echo "========== # of slots : $i ===========" >> time_record.out
#  (time ../../../schedule/simulator -n 150 -H 8 -q 32 -p 0 -b 180.0 -t 0.1 -m   ../../mapFile/mapFile_uni_150_r150/mapFile_uniR150_N150_1.txt -f 0.95 -c 0.477 -A MinPower -T  1.0 -F Kmeans -N $i  -S solution_N150.out -P power_N150.out) 2>> time_record_N150.out
#done

#for (( i = 1; i < 6; i++ )); do
#  echo "========== # of slots : $i ==========="
#  echo "========== # of slots : $i ===========" >> time_record.out
#  (time ../../../schedule/simulator -n 100 -H 5 -q 32 -p 0 -b 180.0 -t 0.1 -m   ../../mapFile/mapFile_uni_100_r150/mapFile_uniR150_N100_1.txt -f 0.95 -c 0.477 -A MinPower -T  1.0 -F Kmeans -N $i  -S solution_N100.out -P power_N100.out) 2>> time_record_N100.out
#done

Fidelity="0.6 0.65 0.7 0.75 0.8"
NumNode=50
NumHead=13
Correlation="0.477"

for item in $Fidelity ; do
  for (( i = 1; i <= 50; i++ )); do
    echo "=========== Fidelity : $item, Index $i ========="
    ../../../schedule/simulator -n $NumNode -H $NumHead -q 32 -p 0 -b 180.0 -t 0.5 -I 10.0 -m  \
      ../../mapFile/mapFile_uni_50_r150/mapFile_uniR150_N50_${i}.txt -f 0.9 -c ${Correlation} -A \
      MinPower  -T 1.0 -F CSFile -N 3 -P data/power_N${NumNode}_H${NumHead}_${i}_F${item}.out -S \
      data/Solution_N${NumNode}_H${NumHead}_${i}_F${item}.out -s cs/data/CS_N50_H13_${i}_F${item}.out
  done
done
