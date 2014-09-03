Params="0.25 0.3 0.35 0.4 0.483"

for ((topoIdx=1;topoIdx<=50;topoIdx=topoIdx+1))
do
  for num in $Params;
  do
    ../../DirectAcessPowerMax_MC/DirectAcessPowerMax_MC 0 180 10.2 ../mapFile/mapFile_uni_195_r500/mapFile_uniR500_N195_${topoIdx}.txt 0.95 ${num} >>  DirectAcessPowerMax195_R500_2_MC.txt
  done
done

