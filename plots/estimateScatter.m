clear all; close all;

fJoin = ...
    dlmread('../runSimulation/runTest/operatorTest/ULSA4b_EstimateJoinGainHN19.txt');

fIsolate = ...
    dlmread('../runSimulation/runTest/operatorTest/ULSA4b_EstimateIsolateGainHN25.txt');
    

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
%idxs = find(vecJoinEstiGain > -500);
figure;
scatter(vecJoinEstiGain,vecJoinRealGain,'*');
title('Joining Operation');
xlabel('Estimated Saved Time (ms)');
ylabel('Real Saved Time (ms)');
%axis([-3 3 -4 1 ])

grid on;

vecIsolateEstiGain = [];
vecIsolateRealGain = []; 

for i = 1:size(fIsolate ,1)
    if fIsolate(i,1) == 1
        numHead = fIsolate(i,3);
        vecIsolateEstiGain = [vecIsolateEstiGain fIsolate(i,2)];
    else
        numHead = fIsolate(i,3);
        vecIsolateRealGain = [vecIsolateRealGain fIsolate(i,2)];
    end
end
%idxs = find(vecJoinEstiGain > -500);
figure;
scatter(vecIsolateEstiGain,vecIsolateRealGain);
title('Isolation Operation');
xlabel('Estimated Saved Time (ms)');
ylabel('Real Saved Time (ms)');
grid on;
