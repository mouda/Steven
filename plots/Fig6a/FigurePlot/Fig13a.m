clear all;
noclu= 26.6148;
MC_colorMap=hsv(16);
DC_colorMap=jet(16);
i=1
hFig = figure(1);
set(hFig, 'Position', [1 1 600 450 ]);
[ totalNum totalInfo MapCR FidRatio  QuantiBits Density BestRound gatherInfo supNum firstEnergy...
    SecondEnergy firstResors secondResors headx] = parseResorsToMatrixByHead...
    ('data/ULSA2i3_All_N195_BW180.0PW0.001_FR0.95_r500.0.txt-1012-multiTopo-20000iters-50times');
mincolsize=returnColNonZeroSize(firstResors);
secondResors2i=secondResors([1:1:mincolsize],[1:1:size(secondResors,2)]);%Compensate
firstResors2i=firstResors([1:1:mincolsize],[1:1:size(firstResors,2)]);
%gatherInfo3i=gatherInfo([1:1:mincolsize],[1:1:size(gatherInfo,2)]);
headx2i=headx;
a=mean(secondResors2i)

[ totalNum totalInfo MapCR FidRatio  QuantiBits Density BestRound gatherInfo supNum firstEnergy...
    SecondEnergy firstResors secondResors headx] = parseResorsToMatrixByHead...
    ('data/ULSA4b2_All_N195_BW180.0PW0.001_FR0.95_r500.0.txt-1012-multiTopo-20000iters-50times');
mincolsize=returnColNonZeroSize(firstResors);
headx4b=headx;
secondResors4b=secondResors([1:1:mincolsize],[1:1:size(secondResors,2)]);
firstResors4b=firstResors([1:1:mincolsize],[1:1:size(firstResors,2)]);
gatherInfo4b=gatherInfo([1:1:mincolsize],[1:1:size(gatherInfo,2)]);


[ totalNum totalInfo MapCR FidRatio  QuantiBits Density BestRound gatherInfo supNum firstEnergy...
    SecondEnergy firstResors secondResors headx] = parseResorsToMatrixByHead...
    ('data/ULSAkmeans4b2_All_N195_BW180.0PW0.001_FR0.95.txt');
mincolsize=returnColNonZeroSize(firstResors);
headxKM4b=headx;
secondResorsKM4b=secondResors([1:1:mincolsize],[1:1:size(secondResors,2)]);
firstResorsKM4b=firstResors([1:1:mincolsize],[1:1:size(firstResors,2)]);
gatherInfoKM4b=gatherInfo([1:1:mincolsize],[1:1:size(gatherInfo,2)]);

[ totalNum totalInfo MapCR FidRatio  QuantiBits Density BestRound gatherInfo supNum firstEnergy...
    SecondEnergy firstResors secondResors headx] = parseResorsToMatrixByHead...
    ('data/ULSAkmeans2i_All_N195_FR0.95.txt');
mincolsize=returnColNonZeroSize(firstResors);
headxKM2i=headx;
secondResorsKM2i=secondResors([1:1:mincolsize],[1:1:size(secondResors,2)]);
firstResorsKM2i=firstResors([1:1:mincolsize],[1:1:size(firstResors,2)]);
gatherInfoKM2i=gatherInfo([1:1:mincolsize],[1:1:size(gatherInfo,2)]);

kk=9;
DDC=16.4416.*ones(1, kk)./noclu;
KM2i=mean(secondResorsKM2i+firstResorsKM2i)./noclu;
KM4b=mean(secondResorsKM4b+firstResorsKM4b)./noclu;
l2i=mean(secondResors2i+firstResors2i)./noclu;
l4b=mean(secondResors4b+firstResors4b)./noclu;
 
 %plot(headxKM2i(1:kk),DDC,'xc--','LineWidth',1.5,'DisplayName','Direct Access-DC','MarkerSize',10);hold on;
 %plot(headxKM2i(1:kk-1),KM2i(1:kk-1),'o--','LineWidth',1.5,'DisplayName','MC Kmeans','Color',MC_colorMap(i+2,:),'MarkerSize',10);hold on;
 %plot(headxKM4b(1:kk-1),KM4b(1:kk-1),'^--','LineWidth',1.5,'DisplayName','DC Kmeans','Color',DC_colorMap(i+2,:),'MarkerSize',10);hold on;
 plot(headx2i(1:kk-1),l2i(1:kk-1),'o-.','LineWidth',1.5,'DisplayName','Two-Tier MC SA','Color',MC_colorMap(i,:),'MarkerSize',10);hold on;
 plot(headx4b(1:kk-1),l4b(1:kk-1),'^-','LineWidth',1.5,'DisplayName','Two-Tier DC SA','Color',DC_colorMap(i,:),'MarkerSize',10);hold on;
% 
% errorbar(headxKM2i,mean(secondResorsKM2i+firstResorsKM2i)./noclu,std((secondResorsKM2i+firstResorsKM2i)./noclu),'o--','LineWidth',1.5,'DisplayName','Kmeans-MC','Color',MC_colorMap(i+2,:),'MarkerSize',10);hold on;
% errorbar(headxKM4b,mean(secondResorsKM4b+firstResorsKM4b)./noclu,std((secondResorsKM4b+firstResorsKM4b)./noclu),'^--','LineWidth',1.5,'DisplayName','Kmeans-DC','Color',DC_colorMap(i+2,:),'MarkerSize',10);hold on;
% errorbar(headx2i,mean(secondResors2i+firstResors2i)./noclu,std((secondResors2i+firstResors2i)./noclu),'o-.','LineWidth',1.5,'DisplayName','Two-Tier MC','Color',MC_colorMap(i,:),'MarkerSize',10);hold on;
% errorbar(headx4b,mean(secondResors4b+firstResors4b)./noclu,std((secondResors4b+firstResors4b)./noclu),'^-','LineWidth',1.5,'DisplayName','Two-Tier DC','Color',DC_colorMap(i,:),'MarkerSize',10);hold on;


title('195 Machines;180 Khz;\lambda=0.95;\eta=0.48');
ylabel('Resource Usage Ratio (R)');
xlabel({'Maximum Allowed Head Number (m)','Two-Tier Data Gathering'});
legend('show');
axis([2 18 0.2 0.77]);
grid on;
