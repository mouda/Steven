#!/bin/bash

#for (( i = 1; i < 6; i++ )); do
#  echo "========== # of slots : $i ==========="
#  echo "========== # of slots : $i ===========" >> time_record.out
#  (time ../../../schedule/simulator -n 50 -H 5 -q 32 -p 0 -b 180.0 -t 0.1 -m   ../../mapFile/mapFile_uni_50_r150/mapFile_uniR150_N50_1.txt -f 0.95 -c 0.477 -A MinPower -T  1.0 -F Kmeans -N $i  -E testEntropy.out -S test.out) 2>> time_record.out
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

for (( i = 1; i < 6; i++ )); do
  echo "========== # of slots : $i ==========="
  echo "========== # of slots : $i ===========" >> time_record.out
  (time ../../../schedule/simulator -n 100 -H 5 -q 32 -p 0 -b 180.0 -t 0.1 -m   ../../mapFile/mapFile_uni_100_r150/mapFile_uniR150_N100_1.txt -f 0.95 -c 0.477 -A MinPower -T  1.0 -F Kmeans -N $i  -S solution_N100.out -P power_N100.out) 2>> time_record_N100.out
done
