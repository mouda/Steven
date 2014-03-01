clear;
close all;
x = 0.0:0.1:0.5;
baselineEntropy = dlmread('entropy_Baseline.out');
branchboundEntropy = dlmread('entropy_Branchbound.out');

str=sprintf('Baseline');
plot(x,mean(baselineEntropy),'*-','LineWidth',2.0,'Color',[0.5 0.5 1],'DisplayName',str,'MarkerSize',10); hold on;
str=sprintf('Branchbound');
plot(x,mean(branchboundEntropy),'^--','LineWidth',2.0,'Color','k','DisplayName',str,'MarkerSize',10); 
legend('show');
xlabel('Time (ms)');
ylabel('Entropy (bits)');
grid on;