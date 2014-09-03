clear; close all;
cellStrEpslion = {'0.001','0.01','0.1','1','10','20'};
numEpslion = [0.001 0.01, 0.1, 1, 10, 20];
bbRation = zeros(1,length(cellStrEpslion));

% Read Branch and Bound Algorithm
for i = 1:length(cellStrEpslion)
    counter = 0;
    entropy = 0;
    for j = 1:50
        
        strFileEpslion = sprintf('data/GE_MaxEntropy_Q8_N50_EPS%s_SC.47_TC.%d.out',cellStrEpslion{i},j);
        strFileBF = sprintf('data_B/GE_BruteForce_Q8_N50_SC.47_TC.%d.out',j);
        s1 = dir(strFileEpslion);
        s2 = dir(strFileBF);
        if exist(strFileEpslion,'file') && s1.bytes ~= 0 && exist(strFileBF,'file') && s2.bytes ~= 0 && dlmread(strFileBF) ~= 0
            data = dlmread(strFileEpslion)/dlmread(strFileBF);
            counter = counter + 1;
            entropy = entropy + data;
        end
    end
    bbRation(i) = entropy/counter;
end

gapRatio = 1 - bbRation;

plot(numEpslion,gapRatio,'-','LineWidth',2,'Color','b');
title('Correlation level $\lambda = 0.477$, $M=50$','interpreter','latex');
ylabel('Performance gap (%)');
xlim([-0.2 20]);
xlabel('Approximation parameter $(\epsilon)$','interpreter','latex');
grid on;