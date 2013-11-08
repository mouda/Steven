clear; close all;

iteration = 1:2000;
vecGuidePayoff = dlmread('data/guideIterObserve.out');
vecNoGuidePayoff = dlmread('data/noGuideIterObserve.out');
DC_colorMap=jet(16);
plot(iteration,vecGuidePayoff(1:2000),'LineWidth',1.5,'DisplayName','Guide Search','Color',DC_colorMap(1,:),'MarkerSize',10); hold on;
plot(iteration,vecNoGuidePayoff(1:2000),'LineWidth',1.5,'DisplayName','No Guide Search','Color',DC_colorMap(3,:),'MarkerSize',10);
legend('show');
grid on;