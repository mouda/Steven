clear; close all;
aryDensity = [ 50 100 150 195];
cellAlg = {'MaxEntropy', 'SumRate','GreedyPhysical','MaxSNR'};
matEntropy = zeros(length(cellAlg),length(aryDensity));

for i = 1:length(cellAlg)
    for j = 1:length(aryDensity)
        counter = 0;
        entropy = 0;
        for k = 1:50
            strFileString = sprintf('data/GE_%s_Q8_N%d_SC.47_TC.1_Idx%d.out',cellAlg{i},aryDensity(j),k);
            s = dir(strFileString);
            if exist(strFileString,'file') && s.bytes ~= 0
                data = dlmread(strFileString);
                entropy = entropy + data;
                counter = counter + 1;
            end
        end
        matEntropy(i,j) = entropy/counter;
    end
end

strMaxEntropy=sprintf('Proposed algorithm');
plot(aryDensity,matEntropy(1,:),'*-','MarkerSize',10, 'LineWidth',2,'Color','r','DisplayName',strMaxEntropy); hold on;
strMaxSumRate=sprintf('Max sum rate');
plot(aryDensity,matEntropy(2,:),'o--','MarkerSize',5,'LineWidth',2,'Color','b','DisplayName',strMaxSumRate); hold on;
strGreedyPhysical=sprintf('GreedyPhysical');
plot(aryDensity,matEntropy(3,:),'x-.','MarkerSize',10,'LineWidth',2,'Color','g','DisplayName',strGreedyPhysical);
strMaxSNR=sprintf('MaxSNR');
plot(aryDensity,matEntropy(4,:),'+:','MarkerSize',8,'LineWidth',2,'Color','k','DisplayName',strMaxSNR);
%ylim([80 400]);
title('Correlation level $\lambda = 0.477$','interpreter','latex');
ylabel('Gathered Entropy (bits)');
xlabel('Number of machines');
legend('show');
grid on;
hold off;