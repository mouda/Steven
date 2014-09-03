clear;
close all;
x = 0.0:0.1:0.5;

baselineTotalEntropy = dlmread('data/TE_Baseline_Q16_N195_SC.47_TC.5.out');
branchboundTotalEntropy = dlmread('data/TE_Branchbound_Q16_N195_SC.47_TC.5.out');
greedyPhysicalTotalEntropy = dlmread('data/TE_GreedyPhysical_Q16_N195_SC.47_TC.5.out');

baselineEntropy = dlmread('data/GE_Baseline_Q16_N195_SC.47_TC.5.out');
branchboundEntropy = dlmread('data/GE_Branchbound_Q16_N195_SC.47_TC.5.out');
greedyPhysical = dlmread('data/GE_GreedyPhysical_Q16_N195_SC.47_TC.5.out');

str=sprintf('Branch and bound algorithm');
plot(x,mean(branchboundEntropy,1),'^--','LineWidth',2.0,'Color','b','DisplayName',str,'MarkerSize',10); hold on;
str=sprintf('GreedyPhysical algorithm');
plot(x,mean(greedyPhysical,1),'^--','LineWidth',2.0,'Color','g','DisplayName',str,'MarkerSize',10);hold on;
str=sprintf('MaxSNR algorithm');
plot(x,mean(baselineEntropy,1),'*-','LineWidth',2.0,'Color','r','DisplayName',str,'MarkerSize',10); 
legend('show');
xlabel('Time (ms)');
ylabel('Entropy (bits)');
grid on;