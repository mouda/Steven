clear;
close all;
x = 0.0:0.1:0.5;
baselineEntropy = dlmread('data/entropy_Baseline_Q_16_N195.out');
branchboundEntropy = dlmread('data/entropy_Branchbound_Q_16_N195.out');
greedyPhysical = dlmread('data/entropy_GreedyPhysical_Q_16_N195.out');

str=sprintf('Branch and bound algorithm');
plot(x,mean(branchboundEntropy,1),'^--','LineWidth',2.0,'Color','b','DisplayName',str,'MarkerSize',10); hold on;
str=sprintf('GreedyPhysical algorithm');
plot(x,mean(greedyPhysical,1),'^--','LineWidth',2.0,'Color','g','DisplayName',str,'MarkerSize',10);hold on;
str=sprintf('Baseline algorithm');
plot(x,mean(baselineEntropy,1),'*-','LineWidth',2.0,'Color','r','DisplayName',str,'MarkerSize',10); 
legend('show');
xlabel('Time (ms)');
ylabel('Entropy (bits)');
grid on;