clear; close all;

numOfIteration = 1000;
iteration = 1:numOfIteration;
vecGuidePayoff = dlmread('data/guideIterObserve.out');
vecNoGuidePayoff = dlmread('data/noGuideIterObserve.out');
DC_colorMap=jet(16);
plot(iteration,vecNoGuidePayoff(1:numOfIteration),'LineWidth',3,'DisplayName','No Guided Local Search','Color','r','MarkerSize',10); hold on;
plot(iteration,vecGuidePayoff(1:numOfIteration),'LineWidth',1.5,'DisplayName','Guided Local Search','Color','b','MarkerSize',10);
ylabel('Payoff (ms)');
xlabel('Iterations');
legend('show');
grid on;