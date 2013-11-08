close all; clear all;
noclu=[2.66148e-005 2.66148e-005 2.66148e-005 2.66148e-005 2.66148e-005];
hFig=figure(1);
set(hFig, 'Position', [1 1 600 450 ]);
CRx= [0.25 0.3 0.35 0.4 0.48 ];
 [ totalNum totalInfo MapCR FidRatio  QuantiBits Density BestRound gatherInfo supNum firstEnergy SecondEnergy firstResors secondResors CRx] ...
   = parseResorsToMatrixByCR_assignedHead( 'data/ULSA4b2_All_N195_BW180.0PW0.001_FR0.95_r500.0.txt',15);
mincolsize=returnColNonZeroSize(firstResors);
Info4b=gatherInfo([1:1:mincolsize],[1:1:size(gatherInfo,2)]);
firR4b=firstEnergy([1:1:mincolsize],[1:1:size(firstResors,2)]);
secR4b=SecondEnergy([1:1:mincolsize],[1:1:size(secondResors,2)]);
totalR4b=firR4b+secR4b;
for j=1:size(totalR4b,2)
FinalRatio4b(:,j)=totalR4b(:,j)./noclu(j);
Final1stRatio4b(:,j)=firR4b(:,j)./noclu(j);
Final2ndRatio4b(:,j)=secR4b(:,j)./noclu(j);
end

[ totalNum totalInfo MapCR FidRatio  QuantiBits Density BestRound gatherInfo supNum firstEnergy SecondEnergy firstResors secondResors CRx]...
    = parseResorsToMatrixByCR_assignedHead( 'data/ULSA2i2_All_N195_BW180.0PW0.001_FR0.95_r500.0.txt',15);
mincolsize=returnColNonZeroSize(firstResors);
secondResors2i=SecondEnergy([1:1:mincolsize],[1:1:size(secondResors,2)])
firstResors2i=firstEnergy([1:1:mincolsize],[1:1:size(secondResors,2)])
gatherInfo2i=gatherInfo([1:1:mincolsize],[1:1:size(gatherInfo,2)])
totalR2i=firstResors2i+secondResors2i;
for j=1:size(totalR2i,2)
    FinalRatio2i(:,j)=totalR2i(:,j)./noclu(j);
    Final1stRatio2i(:,j)=firstResors2i(:,j)./noclu(j);
    Final2ndRatio2i(:,j)=secondResors2i(:,j)./noclu(j);
end

[ totalNum totalInfo MapCR FidRatio  QuantiBits Density BestRound gatherInfo supNum firstEnergy SecondEnergy firstResors secondResors CRx]...
    = parseResorsToMatrixByCR_assignedHead( 'data/UULSAkmeans_DC_All_N195_BW180.0PW0.001_FR0.95_r500.0.txt',11);
mincolsize=returnColNonZeroSize(firstResors);
secondResorsKM4b=SecondEnergy([1:1:mincolsize],[1:1:size(secondResors,2)])
firstResorsKM4b=firstEnergy([1:1:mincolsize],[1:1:size(secondResors,2)])
gatherInfoKM4b=gatherInfo([1:1:mincolsize],[1:1:size(gatherInfo,2)])
totalRKM4b=firstResorsKM4b+secondResorsKM4b;
for j=1:size(totalR2i,2)
    FinalRatioKM4b(:,j)=totalRKM4b(:,j)./noclu(j);
    Final1stRatioKM4b(:,j)=firstResorsKM4b(:,j)./noclu(j);
    Final2ndRatioKM4b(:,j)=secondResorsKM4b(:,j)./noclu(j);
end

[ totalNum totalInfo MapCR FidRatio  QuantiBits Density BestRound gatherInfo supNum firstEnergy SecondEnergy firstResors secondResors CRx]...
    = parseResorsToMatrixByCR_assignedHead( 'data/ULSAkmeans_MC_All_N195_BW180.0PW0.001_FR0.95_r500.0.txt',11);
mincolsize=returnColNonZeroSize(firstResors);
secondResorsKM2i=SecondEnergy([1:1:mincolsize],[1:1:size(secondResors,2)])
firstResorsKM2i=firstEnergy([1:1:mincolsize],[1:1:size(secondResors,2)])
gatherInfoKM2i=gatherInfo([1:1:mincolsize],[1:1:size(gatherInfo,2)])
totalRKM2i=firstResorsKM2i+secondResorsKM2i;
for j=1:size(totalR2i,2)
    FinalRatioKM2i(:,j)=totalRKM2i(:,j)./noclu(j);
    Final1stRatioKM2i(:,j)=firstResorsKM2i(:,j)./noclu(j);
    Final2ndRatioKM2i(:,j)=secondResorsKM2i(:,j)./noclu(j);
end

target(:,1)=mean(FinalRatio4b);
%target(:,2)=mean(FinalRatio2i);
target(:,2)=mean(FinalRatioKM4b);
target(:,3)=mean(FinalRatioKM2i);
index=[1 2 3 4 5];
Obar = bar(index,target,'hist');
set(Obar(1),'FaceColor',[1 0.35 0.35]);
%set(Obar(2),'FaceColor',[0 0 0.7]);
set(Obar(2),'FaceColor',[0.75 1 0.79]);
set(Obar(3),'FaceColor',[1 1 0.9]);
title('195 Machines;180 Khz;\lambda=0.95');
xlabel({'Compression Ratio (\eta)','Two-Tier Data Gathering'});
ylabel('Energy Consumption Ratio');
legend('Data-Centric SA','Data-Centric K-means','Machine-Centric K-means');

set(gca,'XTick',index);
xlabel({'Compression Ratio (\eta)','Two-Tier Data Gathering'});  
set(gca,'XTickLabel',{'0.24';'0.28';'0.40';'0.44';'0.48'});
grid on;
%axis([0.5 1 0 1.2]);
