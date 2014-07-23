clear; close all;

% total time: 100 second
% tier1 time: 90  second
% tier2 time: 10 second
directAccessImagePower = 0.5809;
clusterTier1Power = [0.0855807 9.628064e-02 1.100387e-01 1.283840e-01 1.540706e-01];
clusterTier2Power = [0.00936631 0.000495054 0.000399608 0.00036397  0.000347194];
totalNodes=24;
N = 5;
times = 100;
totalTime = 100000.0;
%read baseline file

% tier1
aryStrTier1Time = {'90000.0' '80000.0' '70000.0' '60000.0' '50000.0'};
tier1Time = [90000.0 80000.0 70000.0 60000.0 50000.0];
baselineClusterTier1Power = zeros(1,length(aryStrTier1Time));
for i=1:length(aryStrTier1Time)
    sumPower = 0;
    countPower = 0;
    for j=1:times
    str = sprintf('data/CS_baseline_%s_%d.out',aryStrTier1Time{i},j);
    [ maxChNum powerMax C2WT payoff SAFac TimeMs1st TimeMs2nd EnergyJ1st EnergyJ2nd clusterStru headList TxPower  ] = parse_Stru_v2 ...
    (str,totalNodes);
        sumPower = sumPower + payoff;
        countPower = countPower + 1;
    end
    baselineClusterTier1Power(i) = sumPower/countPower;
end

% tier2
baselineClusterTier2Power = zeros(1,length(aryStrTier1Time));
for i=1:length(aryStrTier1Time)
    tier2SlotTime = (totalTime - tier1Time(i))/N;
    strTier2Time = sprintf('%d', int64(tier2SlotTime));
    str = sprintf('data/ImageBaseline_%s.out',strTier2Time);
    tmpAry = dlmread(str);
    baselineClusterTier2Power(i) = mean(tmpAry(tmpAry<5));
end

totalVecPowerRatio = (clusterTier1Power+clusterTier2Power)/directAccessImagePower;
xaxis = [90.0 80.0 70.0 60.0 50.0 ];

subplot(2,1,1);
str_baseline = sprintf('Kmeans clustering algorithm');
plot(xaxis, baselineClusterTier1Power, 'x-.','LineWidth',2,'Color','b','DisplayName',str_baseline ); hold on;
str_CSA = sprintf('Proposed clustering algorithm');
plot(xaxis, clusterTier1Power, 'o-','LineWidth',2,'Color','r','DisplayName',str_CSA ); hold on;

xlabel('Tier 1 resource (s)');
ylabel('Tier 1 Power (Watt)');
title('|S|=24, |H|=4, T+N*\tau = 100(s)');
legend('show');
grid on;

subplot(2,1,2);
str_baseline = sprintf('Kmeans clustering algorithm');
plot(xaxis, baselineClusterTier2Power, '*-.','LineWidth',2,'Color', 'b', 'DisplayName',str_baseline); hold on;
str_CSA = sprintf('Proposed clustering algorithm');
plot(xaxis, clusterTier2Power, 'o-','LineWidth',2,'Color', 'r', 'DisplayName',str_CSA); hold on;

xlabel('Tier 1 resource (s)');
ylabel('Tier 2 Power (Watt)');
%title('|S|=24, |H|=4, T+N*\tau = 100(s)');
legend('show');
grid on;
hold off;

