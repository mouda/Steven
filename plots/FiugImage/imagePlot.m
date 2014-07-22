clear; close all;

% total time: 100 second
% tier1 time: 90  second
% tier2 time: 10 second
directAccessImagePower = 0.5809;
clusterTier1Power = [0.0855807 9.628064e-02 1.100387e-01 1.283840e-01 1.540706e-01];
clusterTier2Power = [0.00936631 0.000495054 0.000399608 0.00036397  0.000347194];

totalVecPowerRatio = (clusterTier1Power+clusterTier2Power)/directAccessImagePower;
xaxis = [90.0 80.0 70.0 60.0 50.0 ];
str_CSA = sprintf('Total Power Ratio');
plot(xaxis, totalVecPowerRatio, 'o-','LineWidth',2,'Color','b','DisplayName',str_CSA ); hold on;
%str_CSA = sprintf('Tier1 Power Saving Ratio');
%plot(xaxis, totalVecPowerRatio, 'o-','LineWidth',2,'Color', 'b', 'DisplayName',str_CSA);
xlabel('Tier 1 resource (s)');
ylabel('Total Power ratio');
title('|S|=24, |H|=4, T+N*\tau = 100(s)');
legend('show');
grid on;
hold off;

