clear;
close all;
algorithmName={'Baseline','Branchbound','GreedyPhysical'};
nodeNum=[ 50 100 150 195];

for m = 1:4
    strFileName=sprintf('data/GE_%s_Q16_N%d_SC.47_TC.5.out',algorithmName{1},nodeNum(m));
    baselineEntropy(m) = mean(sum(dlmread(strFileName),2));
end

for m = 1:4
    strFileName=sprintf('data/GE_%s_Q16_N%d_SC.47_TC.5.out',algorithmName{2},nodeNum(m));
    branchboundEntropy(m) = mean(sum(dlmread(strFileName),2));
end

for m = 1:4
    strFileName=sprintf('data/GE_%s_Q16_N%d_SC.47_TC.5.out',algorithmName{3},nodeNum(m));
    greedyPhysical(m) = mean(sum(dlmread(strFileName),2));
end

xAxisLength = length(baselineEntropy);
xAxis = 1:xAxisLength;
str=sprintf('Proposed algorithm');
plot(nodeNum,branchboundEntropy,'^--','LineWidth',2.0,'Color','b','DisplayName',str,'MarkerSize',10); hold on;
str=sprintf('GreedyPhysical algorithm');
plot(nodeNum,greedyPhysical,'^--','LineWidth',2.0,'Color','g','DisplayName',str,'MarkerSize',10);hold on;
str=sprintf('MaxSNR algorithm');
plot(nodeNum,baselineEntropy,'*-','LineWidth',2.0,'Color','r','DisplayName',str,'MarkerSize',10); 

legend('show');
xlabel('Number of machines');
ylabel('Entropy (bits)');
grid on;
