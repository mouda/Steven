clear all;
noclu=26.6148;
MC_colorMap=hsv(16);
DC_colorMap=jet(16);
i=1
hFig = figure(1);
set(hFig, 'Position', [1 1 600 450 ]);
[ totalNum totalInfo MapCR FidRatio  QuantiBits Density BestRound gatherInfo supNum firstEnergy...
    SecondEnergy firstResors secondResors headx] = parseResorsToMatrixByHead...
    ('data/ULSA3i2_All_N195_BW180.0PW0.001_FR0.95.txt');
mincolsize=returnColNonZeroSize(firstResors);
secondResors3i=secondResors([1:1:mincolsize],[1:1:size(secondResors,2)]);%Compensate
firstResors3i=firstResors([1:1:mincolsize],[1:1:size(firstResors,2)]);
%gatherInfo3i=gatherInfo([1:1:mincolsize],[1:1:size(gatherInfo,2)]);
headx3=headx;
a=mean(secondResors3i);
% plot(headx,mean(secondResors3i)./noclu,'^b-','LineWidth',2.5,'DisplayName','Tier-2 Resource(\tau)');hold on;
% plot(headx,mean(firstResors3i)./noclu,'xr--','LineWidth',2.5,'DisplayName','Tier-1 Resource(T_1)');hold on;
% plot(headx,mean(firstResors3i+secondResors3i)./noclu,'^r--','LineWidth',2.5,'DisplayName','Tier-1 Resource(T_1)');hold on;
%plot(headx,mean(gatherInfo3i./(firstResors3i+secondResors3i)),'xr--','LineWidth',2.5,'DisplayName','Tier-1 Resource(T_1)');hold on;

[ totalNum totalInfo MapCR FidRatio  QuantiBits Density BestRound gatherInfo supNum firstEnergy...
    SecondEnergy firstResors secondResors headx] = parseResorsToMatrixByHead...
    ('data/ULSA4b2_All_N195_BW180.0PW0.001_FR0.95_r500.0.txt-1012-multiTopo-20000iters-50times');
mincolsize=returnColNonZeroSize(firstResors);
secondResors4b=secondResors([1:1:mincolsize],[1:1:size(secondResors,2)]);
firstResors4b=firstResors([1:1:mincolsize],[1:1:size(firstResors,2)]);

gatherInfo4b=gatherInfo([1:1:mincolsize],[1:1:size(gatherInfo,2)]);

% plot(headx,mean(secondResors4b)./noclu,'^g-','LineWidth',2.5,'DisplayName','Tier-2 Resource(\tau)');hold on;
% plot(headx,mean(firstResors4b)./noclu,'xk--','LineWidth',2.5,'DisplayName','Tier-1 Resource(T_1)');hold on;
% plot(headx,mean(firstResors4b+secondResors4b)./noclu,'Ok--','LineWidth',2.5,'DisplayName','Tier-1 Resource(T_1)');hold on;
%plot(headx,mean(gatherInfo4b./(firstResors4b+secondResors4b)),'xk--','LineWidth',2.5,'DisplayName','Tier-1 Resource(T_1)');hold on;


plot(headx3,mean(secondResors3i+firstResors3i)./noclu,'*:','LineWidth',1.5,'DisplayName','Two Tiers: Head=m','Color',DC_colorMap(i+4,:),'MarkerSize',10);hold on;
plot(headx,mean(secondResors4b+firstResors4b)./noclu,'*-','LineWidth',1.5,'DisplayName','Two Tiers: Head\leqm','Color',MC_colorMap(i+2,:),'MarkerSize',10);hold on;
plot(headx3,mean(firstResors3i)./noclu,'x:','LineWidth',1.5,'DisplayName','Tier-1: Head=m','Color',DC_colorMap(i,:),'MarkerSize',10);hold on;
plot(headx,mean(firstResors4b)./noclu,'x-','LineWidth',1.5,'DisplayName','Tier-1: Head\leqm','Color',MC_colorMap(i,:),'MarkerSize',10);hold on;
plot(headx3,mean(secondResors3i)./noclu,'d:','LineWidth',1.5,'DisplayName','Tier-2: Head=m','Color',DC_colorMap(i+2,:),'MarkerSize',10);
plot(headx,mean(secondResors4b)./noclu,'d-','LineWidth',1.5,'DisplayName','Tier-2: Head\leqm','Color',MC_colorMap(i+1,:),'MarkerSize',10);


title('195 Machines;180 Khz;\lambda=0.95;\eta=0.48');
ylabel('Resource Usage Ratio (R)');
xlabel({'Maximum Allowed Head Number (m)','Two-Tier Data Gathering'});
legend('show');
axis([3 23 0.0 0.55]);
grid on;
% 

% 
% [ totalNum totalInfo MapCR FidRatio  QuantiBits Density BestRound gatherInfo supNum firstEnergy...
%     SecondEnergy firstResors secondResors headx] = parseResorsToMatrixByHead...
%     ('data/ULSA3gx_DC_All_N195PW0BW180.0_FR0.95.txt');
% mincolsize=returnColNonZeroSize(firstResors);
% secondResors3gx=secondResors([1:1:mincolsize],[1:1:size(secondResors,2)]);
% firstResors3gx=firstResors([1:1:mincolsize],[1:1:size(firstResors,2)]);
% 
% plot(headx,mean(secondResors3gx),'^r-');