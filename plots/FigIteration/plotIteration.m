clear; close all;

numOfIteration = 1000;
iteration = 1:numOfIteration;
vecGuidePayoff = dlmread('data/guideIterObserve.out');
vecNoGuidePayoff = dlmread('data/noGuideIterObserve.out');
DC_colorMap=jet(16);

plot(iteration,vecNoGuidePayoff(1:numOfIteration),'LineWidth',3,'DisplayName','No Guided Local Search','Color','r','MarkerSize',10); hold on;
plot(iteration,vecGuidePayoff(1:numOfIteration),'LineWidth',1.5,'DisplayName','Guided Local Search','Color','b','MarkerSize',10);
title('195 Machines;180 Khz;\lambda=0.95;\eta=0.48');
ylabel('Objective Function (Resource Usage)');
xlabel('Iteration');
legend('show');
grid on;