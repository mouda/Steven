clear all;

hn=[15 11];
noclu=[26.6148 26.6148 26.6148 26.6148 26.6148];
DDC=[ 23.6878 22.9276 21.8165  19.3499 16.4416];
% noclu=[ 1 1 1 1 1];

% directAccess 
rawData = dlmread('data/DirectAcessPowerMax195_R500_2_MC.txt');
x = 1;
for i = 1:50
   noclu(i,:) = rawData(x:x+4,9);
   x = x + 5;
end

MC_colorMap=hsv(16);
DC_colorMap=jet(16);
%hFig=figure(1);
%set(hFig, 'Position', [1 1 600 450 ]);

[ totalNum totalInfo MapCR FidRatio  QuantiBits Density BestRound gatherInfo supNum firstEnergy SecondEnergy firstResors secondResors CRx] ...
   = parseResorsToMatrixByCR_assignedHead( 'data/ULSA4b7_All_N195_BW180.0PW0.001_FR0.95_r500.0.txt',11);
mincolsize=returnColNonZeroSize(firstResors);
Info4b=gatherInfo([1:1:mincolsize],[1:1:size(gatherInfo,2)]);
firR4b=firstResors([1:1:mincolsize],[1:1:size(firstResors,2)]);
secR4b=secondResors([1:1:mincolsize],[1:1:size(secondResors,2)]);
totalR4b=firR4b+secR4b;
for j=1:size(totalR4b,2)
FinalRatio4b(:,j)=totalR4b(:,j)./noclu(:,j);
Final1stRatio4b(:,j)=firR4b(:,j)./noclu(:,j);
Final2ndRatio4b(:,j)=secR4b(:,j)./noclu(:,j);
end
% errorbar(CRx,mean(Final1stRatio3i),std(Final1stRatio3i),'^--','LineWidth',2.5,'Color',DC_colorMap(i+1,:),'DisplayName','Data-Centric: Tier-1'); hold on;
% errorbar(CRx,mean(Final2ndRatio3i),std(Final2ndRatio3i),'^-.','LineWidth',2.5,'Color',DC_colorMap(i+2,:),'DisplayName','Data-Centric: Tier-2'); hold on;

%errorbar(CRx,mean(Final1stRatio3i),std(Final1stRatio3i),'^--','LineWidth',2.5,'Color',DC_colorMap(i+1,:),'DisplayName','DC:1st'); hold on;
%errorbar(CRx,mean(Final2ndRatio3i),std(Final2ndRatio3i),'^-.','LineWidth',2.5,'Color',DC_colorMap(i+2,:),'DisplayName','DC:2nd'); hold on;


%plot(CRx,mean(totalR3i),'^-','LineWidth',2.5,'Color',DC_colorMap(i,:),'DisplayName',str); hold on;

[ totalNum totalInfo MapCR FidRatio  QuantiBits Density BestRound gatherInfo supNum firstEnergy SecondEnergy firstResors secondResors CRx]...
    = parseResorsToMatrixByCR_assignedHead( 'data/ULSAkmeans_DC_All_N195_BW180.0PW0.001_FR0.95_r500.0.txt',hn(2));
mincolsize=returnColNonZeroSize(firstResors);
secondResorsKM4b=secondResors([1:1:mincolsize],[1:1:size(secondResors,2)]);
firstResorsKM4b=firstResors([1:1:mincolsize],[1:1:size(secondResors,2)]);
gatherInfoKM4b=gatherInfo([1:1:mincolsize],[1:1:size(gatherInfo,2)]);
totalRKM4b=firstResorsKM4b+secondResorsKM4b;
for j=1:size(totalRKM4b,2)
    FinalRatioKM4b(:,j)=totalRKM4b(:,j)./noclu(:,j);
    Final1stRatioKM4b(:,j)=firstResorsKM4b(:,j)./noclu(:,j);
    Final2ndRatioKM4b(:,j)=secondResorsKM4b(:,j)./noclu(:,j);
end

[ totalNum totalInfo MapCR FidRatio  QuantiBits Density BestRound gatherInfo supNum firstEnergy SecondEnergy firstResors secondResors CRx]...
    = parseResorsToMatrixByCR_assignedHead( 'data/ULSAkmeans_MC_All_N195_BW180.0PW0.001_FR0.95_r500.0.txt',hn(2));
mincolsize=returnColNonZeroSize(firstResors);
secondResorsKM2i=secondResors([1:1:mincolsize],[1:1:size(secondResors,2)]);
firstResorsKM2i=firstResors([1:1:mincolsize],[1:1:size(secondResors,2)]);
gatherInfoKM2i=gatherInfo([1:1:mincolsize],[1:1:size(gatherInfo,2)]);
totalRKM2i=firstResorsKM2i+secondResorsKM2i;
for j=1:size(totalRKM2i,2)
    FinalRatioKM2i(:,j)=totalRKM2i(:,j)./noclu(:,j);
    Final1stRatioKM2i(:,j)=firstResorsKM2i(:,j)./noclu(:,j);
    Final2ndRatioKM2i(:,j)=secondResorsKM2i(:,j)./noclu(:,j);
end

fig = figure(1);
figH = subplot(2,1,1,'Parent',fig);
str=sprintf('Direct Access-DC');
str=sprintf('MC Kmeans (M=%d)',hn(2));
errorbar(figH, CRx,mean(FinalRatioKM2i),std(FinalRatioKM2i),'o-.','LineWidth',1.5,'Color','r','DisplayName',str,'MarkerSize',10); hold(figH,'on');
str=sprintf('DC Kmeans (M=%d)',hn(2));
errorbar(figH, CRx,mean(FinalRatioKM4b),std(FinalRatioKM4b),'^--','LineWidth',1.5,'Color','k','DisplayName',str,'MarkerSize',10); hold(figH,'on');
str=sprintf('Two-Tier DC SA');
errorbar(figH,CRx,mean(FinalRatio4b),std(FinalRatio4b),'^-','LineWidth',1.5,'Color','b','DisplayName',str,'MarkerSize',10); 



grid on;
legend('show');
ylabel('Resource Usage Ratio (R)');
title('195 Machines;180 Khz;\lambda=0.95');
hold off
ylim([0.15 1.0]);
xlim([0.2 0.5]);
rawData = dlmread('data/DirectAcessPowerMax195_R500_2_MC.txt');
x = 1;
for i = 1:50
   noclu(i,:) = rawData(x:x+4,10);
   x = x + 5;
end

CRx= [0.25 0.3 0.35 0.4 0.48 ];
 [ totalNum totalInfo MapCR FidRatio  QuantiBits Density BestRound gatherInfo supNum firstEnergy SecondEnergy firstResors secondResors CRx] ...
   = parseResorsToMatrixByCR_assignedHead( 'data/ULSA4b7_All_N195_BW180.0PW0.001_FR0.95_r500.0.txt',11);
mincolsize=returnColNonZeroSize(firstResors);
Info4b=gatherInfo([1:1:mincolsize],[1:1:size(gatherInfo,2)]);
firR4b=firstEnergy([1:1:mincolsize],[1:1:size(firstResors,2)]);
secR4b=SecondEnergy([1:1:mincolsize],[1:1:size(secondResors,2)]);
totalR4b=firR4b+secR4b;
for j=1:size(totalR4b,2)
FinalRatio4b(:,j)=totalR4b(:,j)./noclu(:,j);
Final1stRatio4b(:,j)=firR4b(:,j)./noclu(:,j);
Final2ndRatio4b(:,j)=secR4b(:,j)./noclu(:,j);
end

[ totalNum totalInfo MapCR FidRatio  QuantiBits Density BestRound gatherInfo supNum firstEnergy SecondEnergy firstResors secondResors CRx]...
    = parseResorsToMatrixByCR_assignedHead( 'data/ULSAkmeans_DC_All_N195_BW180.0PW0.001_FR0.95_r500.0.txt',11);
mincolsize=returnColNonZeroSize(firstResors);
secondResorsKM4b=SecondEnergy([1:1:mincolsize],[1:1:size(secondResors,2)])
firstResorsKM4b=firstEnergy([1:1:mincolsize],[1:1:size(secondResors,2)])
gatherInfoKM4b=gatherInfo([1:1:mincolsize],[1:1:size(gatherInfo,2)])
totalRKM4b=firstResorsKM4b+secondResorsKM4b;
for j=1:size(totalRKM4b,2)
    FinalRatioKM4b(:,j)=totalRKM4b(:,j)./noclu(:,j);
    Final1stRatioKM4b(:,j)=firstResorsKM4b(:,j)./noclu(:,j);
    Final2ndRatioKM4b(:,j)=secondResorsKM4b(:,j)./noclu(:,j);
end

[ totalNum totalInfo MapCR FidRatio  QuantiBits Density BestRound gatherInfo supNum firstEnergy SecondEnergy firstResors secondResors CRx]...
    = parseResorsToMatrixByCR_assignedHead( 'data/ULSAkmeans_MC_All_N195_BW180.0PW0.001_FR0.95_r500.0.txt',11);
mincolsize=returnColNonZeroSize(firstResors);
secondResorsKM2i=SecondEnergy([1:1:mincolsize],[1:1:size(secondResors,2)])
firstResorsKM2i=firstEnergy([1:1:mincolsize],[1:1:size(secondResors,2)])
gatherInfoKM2i=gatherInfo([1:1:mincolsize],[1:1:size(gatherInfo,2)])
totalRKM2i=firstResorsKM2i+secondResorsKM2i;
for j=1:size(totalRKM2i,2)
    FinalRatioKM2i(:,j)=totalRKM2i(:,j)./noclu(:,j);
    Final1stRatioKM2i(:,j)=firstResorsKM2i(:,j)./noclu(:,j);
    Final2ndRatioKM2i(:,j)=secondResorsKM2i(:,j)./noclu(:,j);
end



figL = subplot(2,1,2,'Parent',fig);

target(:,1)=mean(FinalRatio4b);
target(:,2)=mean(FinalRatioKM4b);
target(:,3)=mean(FinalRatioKM2i);
index=[1 2 3 4 5];
Obar = bar(figL, index,target,'hist');
%axis([0.2 0.5 0.1 1.2]);
ylim([0.1 1.0]);
%xlim([0.2 0.5]);
set(Obar(1),'FaceColor',[1 0.35 0.35]);
set(Obar(2),'FaceColor',[0 1 0.7]);
set(Obar(3),'FaceColor',[1 1 0.9]);

xlabel({'Compression Ratio (\eta)','Two-Tier Data Gathering'});
ylabel('Energy Consumption Ratio');
legend('Two-Tier Data Gathering','DC Kmeans (M=11)','MC Kmeans (M=11)');
set(gca,'XTickLabel',{'0.25';'0.3';'0.35';'0.4';'0.48'});
grid on;
set(gca,'XTick',index);%axis([0.2 0.5 0.1 1.2]);
hold off;


