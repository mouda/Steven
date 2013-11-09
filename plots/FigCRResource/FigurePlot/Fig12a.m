clear all;

hn=[15 11];
noclu=[26.6148 26.6148 26.6148 26.6148 26.6148];
DDC=[ 23.6878 22.9276 21.8165  19.3499 16.4416];
% noclu=[ 1 1 1 1 1];

MC_colorMap=hsv(16);
DC_colorMap=jet(16);
hFig=figure(1);
set(hFig, 'Position', [1 1 600 450 ]);
for i=1:1
 [ totalNum totalInfo MapCR FidRatio  QuantiBits Density BestRound gatherInfo supNum firstEnergy SecondEnergy firstResors secondResors CRx] ...
   = parseResorsToMatrixByCR_assignedHead( 'data/ULSA4b2_All_N195_BW180.0PW0.001_FR0.95_r500.0.txt',15);
mincolsize=returnColNonZeroSize(firstResors);
Info4b=gatherInfo([1:1:mincolsize],[1:1:size(gatherInfo,2)]);
firR4b=firstResors([1:1:mincolsize],[1:1:size(firstResors,2)]);
secR4b=secondResors([1:1:mincolsize],[1:1:size(secondResors,2)]);
totalR4b=firR4b+secR4b;
for j=1:size(totalR4b,2)
FinalRatio4b(:,j)=totalR4b(:,j)./noclu(j);
Final1stRatio4b(:,j)=firR4b(:,j)./noclu(j);
Final2ndRatio4b(:,j)=secR4b(:,j)./noclu(j);
end
% errorbar(CRx,mean(Final1stRatio3i),std(Final1stRatio3i),'^--','LineWidth',2.5,'Color',DC_colorMap(i+1,:),'DisplayName','Data-Centric: Tier-1'); hold on;
% errorbar(CRx,mean(Final2ndRatio3i),std(Final2ndRatio3i),'^-.','LineWidth',2.5,'Color',DC_colorMap(i+2,:),'DisplayName','Data-Centric: Tier-2'); hold on;

%errorbar(CRx,mean(Final1stRatio3i),std(Final1stRatio3i),'^--','LineWidth',2.5,'Color',DC_colorMap(i+1,:),'DisplayName','DC:1st'); hold on;
%errorbar(CRx,mean(Final2ndRatio3i),std(Final2ndRatio3i),'^-.','LineWidth',2.5,'Color',DC_colorMap(i+2,:),'DisplayName','DC:2nd'); hold on;


%plot(CRx,mean(totalR3i),'^-','LineWidth',2.5,'Color',DC_colorMap(i,:),'DisplayName',str); hold on;

[ totalNum totalInfo MapCR FidRatio  QuantiBits Density BestRound gatherInfo supNum firstEnergy SecondEnergy firstResors secondResors CRx]...
    = parseResorsToMatrixByCR_assignedHead( 'data/ULSA2i2_All_N195_BW180.0PW0.001_FR0.95_r500.0.txt',15);
mincolsize=returnColNonZeroSize(firstResors);
secondResors2i=secondResors([1:1:mincolsize],[1:1:size(secondResors,2)]);
firstResors2i=firstResors([1:1:mincolsize],[1:1:size(secondResors,2)]);
gatherInfo2i=gatherInfo([1:1:mincolsize],[1:1:size(gatherInfo,2)]);
totalR2i=firstResors2i+secondResors2i;
for j=1:size(totalR2i,2)
    FinalRatio2i(:,j)=totalR2i(:,j)./noclu(j);
    Final1stRatio2i(:,j)=firstResors2i(:,j)./noclu(j);
    Final2ndRatio2i(:,j)=secondResors2i(:,j)./noclu(j);
end
a=mean(FinalRatio2i);
b=std(FinalRatio2i);

[ totalNum totalInfo MapCR FidRatio  QuantiBits Density BestRound gatherInfo supNum firstEnergy SecondEnergy firstResors secondResors CRx]...
    = parseResorsToMatrixByCR_assignedHead( 'data/ULSAkmeans_DC_All_N195_BW180.0PW0.001_FR0.95_r500.0.txt',hn(2));
mincolsize=returnColNonZeroSize(firstResors);
secondResorsKM4b=secondResors([1:1:mincolsize],[1:1:size(secondResors,2)]);
firstResorsKM4b=firstResors([1:1:mincolsize],[1:1:size(secondResors,2)]);
gatherInfoKM4b=gatherInfo([1:1:mincolsize],[1:1:size(gatherInfo,2)]);
totalRKM4b=firstResorsKM4b+secondResorsKM4b;
for j=1:size(totalR2i,2)
    FinalRatioKM4b(:,j)=totalRKM4b(:,j)./noclu(j);
    Final1stRatioKM4b(:,j)=firstResorsKM4b(:,j)./noclu(j);
    Final2ndRatioKM4b(:,j)=secondResorsKM4b(:,j)./noclu(j);
end

[ totalNum totalInfo MapCR FidRatio  QuantiBits Density BestRound gatherInfo supNum firstEnergy SecondEnergy firstResors secondResors CRx]...
    = parseResorsToMatrixByCR_assignedHead( 'data/ULSAkmeans_MC_All_N195_BW180.0PW0.001_FR0.95_r500.0.txt',hn(2));
mincolsize=returnColNonZeroSize(firstResors);
secondResorsKM2i=secondResors([1:1:mincolsize],[1:1:size(secondResors,2)]);
firstResorsKM2i=firstResors([1:1:mincolsize],[1:1:size(secondResors,2)]);
gatherInfoKM2i=gatherInfo([1:1:mincolsize],[1:1:size(gatherInfo,2)]);
totalRKM2i=firstResorsKM2i+secondResorsKM2i;
for j=1:size(totalR2i,2)
    FinalRatioKM2i(:,j)=totalRKM2i(:,j)./noclu(j);
    Final1stRatioKM2i(:,j)=firstResorsKM2i(:,j)./noclu(j);
    Final2ndRatioKM2i(:,j)=secondResorsKM2i(:,j)./noclu(j);
end


str=sprintf('Direct Access-DC');
%plot(CRx,DDC./noclu,'ro-','LineWidth',2.5,'Color',DC_colorMap(i+4,:),'DisplayName',str,'MarkerSize',10); hold on;
str=sprintf('MC Kmeans (m=%d)',hn(2));
errorbar(CRx,mean(FinalRatioKM2i),std(FinalRatioKM2i),'o--','LineWidth',1.5,'Color',MC_colorMap(i+2,:),'DisplayName',str,'MarkerSize',10); hold on;
str=sprintf('DC Kmeans (m=%d)',hn(2));
errorbar(CRx,mean(FinalRatioKM4b),std(FinalRatioKM4b),'^--','LineWidth',1.5,'Color',DC_colorMap(2+i,:),'DisplayName',str,'MarkerSize',10); hold on;
%str=sprintf('Two-Tier MC SA');
%errorbar(CRx,a,b,'o-','LineWidth',1.5,'Color',MC_colorMap(i,:),'DisplayName',str,'MarkerSize',10); hold on;
str=sprintf('Two-Tier DC SA');
errorbar(CRx,mean(FinalRatio4b),std(FinalRatio4b),'^-','LineWidth',2.5,'Color',DC_colorMap(i,:),'DisplayName',str,'MarkerSize',10); hold on;

%errorbar(CRx,mean(Final1stRatio2i),std(Final1stRatio2i),'o--','LineWidth',2.5,'Color',MC_colorMap(i+1,:),'DisplayName','MC:1st'); hold on;
%errorbar(CRx,mean(Final2ndRatio2i),std(Final2ndRatio2i),'o-.','LineWidth',2.5,'Color',MC_colorMap(i+2,:),'DisplayName','MC:2nd'); hold on;


end
grid on;
legend('show');
xlabel({'Compression Ratio (\eta)','Two-Tier Data Gathering'});
ylabel('Resource Usage Ratio (R)');
title('195 Machines;180 Khz;\lambda=0.95');