clear all; close all;
fIsolate = ...
    dlmread('../runSimulation/runTest/observation/ULSA4b_EstimateIsolateGainHN25.txt');
fJoin = ...
    dlmread('../runSimulation/runTest/observation/ULSA4b_EstimateJoinGainHN25.txt');

maxHeadNum = 25;
matIsoRealGain = zeros(maxHeadNum,size(fIsolate,1));
matIsoEstiGain = zeros(maxHeadNum,size(fIsolate,1));
matJoinRealGain = zeros(maxHeadNum,size(fJoin,1));
matJoinEstiGain = zeros(maxHeadNum, size(fJoin,1));
for i = 1:size(fIsolate,1)
    if fIsolate(i,1) == 1
        numHead = fIsolate(i,3);
        matIsoEstiGain(numHead, find(matIsoEstiGain(numHead,:) == 0,1) ) = fIsolate(i,2);
    else
        matIsoRealGain(numHead, find(matIsoRealGain(numHead,:) == 0,1) ) = fIsolate(i,2);
    end
end

for i = 1:size(fJoin,1)
    if fJoin(i,1) == 0
        numHead = fJoin(i,3);
        matJoinEstiGain(numHead, find(matJoinEstiGain(numHead,:) == 0,1) ) = fJoin(i,2);
    else
        numHead = fJoin(i,3);
        matJoinRealGain(numHead, find(matJoinRealGain(numHead,:) == 0,1) ) = fJoin(i,2);
    end
end

colorMap=jet(4);


meanIsoRealGain = zeros(maxHeadNum,1);
meanIsoEstiGain = zeros(maxHeadNum,1);
meanJoinRealGain = zeros(maxHeadNum,1);
meanJoinEstiGain = zeros(maxHeadNum,1);

stdIsoRealGain = zeros(maxHeadNum,1);
stdIsoEstiGain = zeros(maxHeadNum,1);
stdJoinRealGain = zeros(maxHeadNum,1);
stdJoinEstiGain = zeros(maxHeadNum,1);


for i = 1:maxHeadNum
    meanIsoRealGain(i) = -1*mean(matIsoRealGain(i, matIsoRealGain(i,:) ~= 0));
    meanIsoEstiGain(i) = -1*mean(matIsoEstiGain(i, matIsoEstiGain(i,:) ~= 0));
    meanJoinRealGain(i) = -1*mean(matJoinRealGain(i, matJoinRealGain(i,:) ~= 0));
    meanJoinEstiGain(i) = -1*mean(matJoinEstiGain(i, matJoinEstiGain(i,:) ~= 0));

    stdIsoRealGain(i) = std(matIsoRealGain(i, matIsoRealGain(i,:) ~= 0));
    stdIsoEstiGain(i) = std(matIsoEstiGain(i, matIsoEstiGain(i,:) ~= 0));
    stdJoinRealGain(i) = std(matJoinRealGain(i, matJoinRealGain(i,:) ~= 0));
    stdJoinEstiGain(i) = std(matJoinEstiGain(i, matJoinEstiGain(i,:) ~= 0));
end
figure;
errorbar(1:25,meanIsoRealGain,stdIsoRealGain,'^-','LineWidth',2.5, 'Color', colorMap(1,:), 'DisplayName','Isolate Real'); hold on;
errorbar(1:25,meanIsoEstiGain,stdIsoEstiGain,'^-','LineWidth',2.5, 'Color', colorMap(2,:), 'DisplayName','Isolate Estimate'); hold off;
grid on;

figure;
errorbar(1:25,meanJoinRealGain,stdJoinRealGain,'^-','LineWidth',2.5, 'Color', colorMap(3,:), 'DisplayName','Join Real'); hold on;
errorbar(1:25,meanJoinEstiGain,stdJoinEstiGain,'^-','LineWidth',2.5, 'Color', colorMap(4,:), 'DisplayName','Join Estimate'); hold off;
grid on;