clear;
close all;
x = 0.0:0.1:1.1;
baselineEntropy = dlmread('data/entropy_Baseline_N50.out');
branchboundEntropy = dlmread('data/entropy_Branchbound_N50.out');


str=sprintf('Branch and bound algorithm');
plot(x,mean(branchboundEntropy,1),'^--','LineWidth',2.0,'Color','b','DisplayName',str,'MarkerSize',10); hold on;
str=sprintf('Baseline algorithm');
plot(x,mean(baselineEntropy,1),'*-','LineWidth',2.0,'Color','r','DisplayName',str,'MarkerSize',10); 
legend('show');
xlabel('Time (ms)');
ylabel('Entropy (bits)');
grid on;