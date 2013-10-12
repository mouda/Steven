clear all; close all;
noclu=16.4416+(30.7786-24.0785);

hFig = figure(1);
set(hFig, 'Position', [1 1 600 450 ]);
DC_colorMap=jet(16);

[ totalNum totalInfo MapCR FidRatio  QuantiBits Density BestRound gatherInfo supNum firstEnergy SecondEnergy ... 
    firstResors secondResors headx finalHead] = parseResorsToMatrixByHead_4b('data/ULSA4b2_All_N195_BW180.0PW0.001_FR0.95_r500.0.txt-1012-80000iters-12times');
mincolsize=returnColNonZeroSize(firstResors);
secondResors4b195=secondResors([1:1:mincolsize],[1:1:size(secondResors,2)]);
firstResors4b195=firstResors([1:1:mincolsize],[1:1:size(firstResors,2)]);
gatherInfo4b195=gatherInfo([1:1:mincolsize],[1:1:size(gatherInfo,2)]);
finalHead195=finalHead;
[ totalNum totalInfo MapCR FidRatio  QuantiBits Density BestRound gatherInfo supNum firstEnergy SecondEnergy ... 
    firstResors secondResors headx finalHead] = parseResorsToMatrixByHead_4b('data/ULSA4b2_All_N150_BW180.0PW0.001_FR0.95_r500.0.txt-1012-80000iters-12times');
mincolsize=returnColNonZeroSize(firstResors);
secondResors4b150=secondResors([1:1:mincolsize],[1:1:size(secondResors,2)]);
firstResors4b150=firstResors([1:1:mincolsize],[1:1:size(firstResors,2)]);
gatherInfo4b150=gatherInfo([1:1:mincolsize],[1:1:size(gatherInfo,2)]);
finalHead150=finalHead;

[ totalNum totalInfo MapCR FidRatio  QuantiBits Density BestRound gatherInfo supNum firstEnergy SecondEnergy ... 
    firstResors secondResors headx finalHead] = parseResorsToMatrixByHead_4b('data/ULSA4b2_All_N50_BW180.0PW0.001_FR0.95_r500.0.txt-1012-80000iters-12times');
mincolsize=returnColNonZeroSize(firstResors);
secondResors4b50=secondResors([1:1:mincolsize],[1:1:size(secondResors,2)]);
firstResors4b50=firstResors([1:1:mincolsize],[1:1:size(firstResors,2)]);
gatherInfo4b50=gatherInfo([1:1:mincolsize],[1:1:size(gatherInfo,2)]);
finalHead50=finalHead;

totalRseource4b50 = firstResors4b50+secondResors4b50;
totalResource4b150 = firstResors4b150+secondResors4b150;
totalResource4b195 = firstResors4b195+secondResors4b195;


[result, idx] = min(firstResors4b50+secondResors4b50,[],1);
bestFinalHead50 = zeros(1,size(finalHead50,2));
for i = 1:length(idx) 
    bestFinalHead50(1,i) = finalHead(idx(i),i);
end

[result, idx] = min(firstResors4b150+secondResors4b150,[],1);
bestFinalHead150 = zeros(1,size(finalHead150,2));
for i = 1:length(idx) 
    bestFinalHead150(1,i) = finalHead(idx(i),i);
end

[result, idx] = min(firstResors4b195+secondResors4b195,[],1);
bestFinalHead195 = zeros(1,size(finalHead195,2));
for i = 1:length(idx) 
    bestFinalHead195(1,i) = finalHead(idx(i),i);
end


 errorbar(headx,mean(totalResource4b195),std(totalResource4b195),'*-b','LineWidth',1.5,'DisplayName','Head\leqm','MarkerSize',10);hold on;
 errorbar(headx,mean(totalResource4b150),std(totalResource4b150),'^--r','LineWidth',1.5,'DisplayName','Head\leqm','MarkerSize',10);hold on;
 errorbar(headx,mean(totalRseource4b50),std(totalRseource4b50),'kx-.','LineWidth',1.5,'DisplayName','Head\leqm','MarkerSize',10);hold on;
%plot(headx,bestFinalHead195,'*-b','LineWidth',1.5,'DisplayName','Head\leqm (195 Machines)','MarkerSize',10);hold on;
%plot(headx,bestFinalHead150,'^--r','LineWidth',1.5,'DisplayName','Head\leqm (150 Machines)','MarkerSize',10);hold on;
%plot(headx,bestFinalHead50,'kx-.','LineWidth',1.5,'DisplayName','Head\leqm (50 Machines)','MarkerSize',10);hold on;



%plot(headx,mean(finalHead),'^-.','LineWidth',1.5,'DisplayName','Head\leqm','MarkerSize',10);hold on;

plot (headx,headx,'kx-','LineWidth',1.5,'DisplayName','Head=m');

title('180 Khz;\lambda=0.95;\eta=0.48');
ylabel('Final Converged Head Number ');
xlabel({'Maximum Allowed Head Number (m)','Two-tier Data Gathering'});
legend('show');
grid on;
axis([3 23 4 23 ]);
