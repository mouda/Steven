clear all; close all;

fJoin = ...
    dlmread('../runSimulation/runTest/operatorTest/ULSA4b_EstimateJoinGainHN25.txt');


vecJoinEstiGain = [];
vecJoinRealGain = [];

for i = 1:size(fJoin,1)
    if fJoin(i,1) == 0
        numHead = fJoin(i,3);
        vecJoinEstiGain = [vecJoinEstiGain fJoin(i,2)];
    else
        numHead = fJoin(i,3);
        vecJoinRealGain = [vecJoinRealGain fJoin(i,2)];
    end
end
scatter(vecJoinEstiGain,vecJoinRealGain);
grid on;