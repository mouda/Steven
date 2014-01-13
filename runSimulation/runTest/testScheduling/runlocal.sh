../../../schedule/simulator -n 50 -H 10 -q 8 -p 0 -b 180.0 -t 0.01 -m ../../mapFile/mapFile_uni_50_r500/mapFile_uniR500_N50_1.txt -c 0.477 -f 0.95 -G Branchbound 
../../../schedule/simulator -n 10 -H 2 -q 8 -p 0 -b 180.0 -t 0.01 -m ../../mapFile/mapFile_uni_10_r150/mapFile_uniR150_N10_1.txt -f 0.95 -c 0.477 -G Branchbound
../../../schedule/simulator -n 10 -H 2 -q 8 -p 0 -b 180.0 -t 0.01 -m ../../mapFile/mapFile_uni_10_r150/mapFile_uniR150_N10_1.txt -f 0.95 -c 0.477 -G Baseline
