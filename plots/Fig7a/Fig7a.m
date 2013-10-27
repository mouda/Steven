clear all; close all;
noclu=16.4416+(30.7786-24.0785);

hFig = figure(1);
set(hFig, 'Position', [1 1 600 450 ]);
DC_colorMap=jet(16);

[ totalNum totalInfo MapCR FidRatio  QuantiBits Density BestRound gatherInfo supNum firstEnergy SecondEnergy ... 
    firstResors secondResors headx finalHead] = parseResorsToMatrixByHead_4b('data/ULSA4b2_All_N50_BW180.0PW0.001_FR0.95_r500.0.txt-1012-multiTopo-20000iters-50times');
mincolsize=returnColNonZeroSize(firstResors);
secondResors4b50=secondResors([1:1:mincolsize],[1:1:size(secondResors,2)]);
firstResors4b50=firstResors([1:1:mincolsize],[1:1:size(firstResors,2)]);
gatherInfo4b50=gatherInfo([1:1:mincolsize],[1:1:size(gatherInfo,2)]);
finalHead50=finalHead;

[ totalNum totalInfo MapCR FidRatio  QuantiBits Density BestRound gatherInfo supNum firstEnergy SecondEnergy ... 
    firstResors secondResors headx finalHead] = parseResorsToMatrixByHead_4b('data/ULSA4b2_All_N75_BW180.0PW0.001_FR0.95_r500.0.txt-1012-multiTopo-20000iters-50times');
mincolsize=returnColNonZeroSize(firstResors);
secondResors4b75=secondResors([1:1:mincolsize],[1:1:size(secondResors,2)]);
firstResors4b75=firstResors([1:1:mincolsize],[1:1:size(firstResors,2)]);
gatherInfo4b75=gatherInfo([1:1:mincolsize],[1:1:size(gatherInfo,2)]);
finalHead75=finalHead;

[ totalNum totalInfo MapCR FidRatio  QuantiBits Density BestRound gatherInfo supNum firstEnergy SecondEnergy ... 
    firstResors secondResors headx finalHead] = parseResorsToMatrixByHead_4b('data/ULSA4b2_All_N100_BW180.0PW0.001_FR0.95_r500.0.txt-1012-multiTopo-20000iters-50times');
mincolsize=returnColNonZeroSize(firstResors);
secondResors4b100=secondResors([1:1:mincolsize],[1:1:size(secondResors,2)]);
firstResors4b100=firstResors([1:1:mincolsize],[1:1:size(firstResors,2)]);
gatherInfo4b100=gatherInfo([1:1:mincolsize],[1:1:size(gatherInfo,2)]);
finalHead100=finalHead;

[ totalNum totalInfo MapCR FidRatio  QuantiBits Density BestRound gatherInfo supNum firstEnergy SecondEnergy ... 
    firstResors secondResors headx finalHead] = parseResorsToMatrixByHead_4b('data/ULSA4b2_All_N125_BW180.0PW0.001_FR0.95_r500.0.txt-1012-multiTopo-20000iters-50times');
mincolsize=returnColNonZeroSize(firstResors);
secondResors4b125=secondResors([1:1:mincolsize],[1:1:size(secondResors,2)]);
firstResors4b125=firstResors([1:1:mincolsize],[1:1:size(firstResors,2)]);
gatherInfo4b125=gatherInfo([1:1:mincolsize],[1:1:size(gatherInfo,2)]);
finalHead125=finalHead;

[ totalNum totalInfo MapCR FidRatio  QuantiBits Density BestRound gatherInfo supNum firstEnergy SecondEnergy ... 
    firstResors secondResors headx finalHead] = parseResorsToMatrixByHead_4b('data/ULSA4b2_All_N150_BW180.0PW0.001_FR0.95_r500.0.txt-1012-multiTopo-20000iters-50times');
mincolsize=returnColNonZeroSize(firstResors);
secondResors4b150=secondResors([1:1:mincolsize],[1:1:size(secondResors,2)]);
firstResors4b150=firstResors([1:1:mincolsize],[1:1:size(firstResors,2)]);
gatherInfo4b150=gatherInfo([1:1:mincolsize],[1:1:size(gatherInfo,2)]);
finalHead150=finalHead;

[ totalNum totalInfo MapCR FidRatio  QuantiBits Density BestRound gatherInfo supNum firstEnergy SecondEnergy ... 
    firstResors secondResors headx finalHead] = parseResorsToMatrixByHead_4b('data/ULSA4b2_All_N195_BW180.0PW0.001_FR0.95_r500.0.txt-1012-multiTopo-20000iters-50times');
mincolsize=returnColNonZeroSize(firstResors);
secondResors4b195=secondResors([1:1:mincolsize],[1:1:size(secondResors,2)]);
firstResors4b195=firstResors([1:1:mincolsize],[1:1:size(firstResors,2)]);
gatherInfo4b195=gatherInfo([1:1:mincolsize],[1:1:size(gatherInfo,2)]);
finalHead195=finalHead;



[ totalNum totalInfo MapCR FidRatio  QuantiBits Density BestRound gatherInfo supNum firstEnergy SecondEnergy ... 
    firstResors secondResors headx finalHead] = parseResorsToMatrixByHead_4b('data/ULSA4b4_All_N50_BW180.0PW0.001_FR0.95_r500.0.txt');
mincolsize=returnColNonZeroSize(firstResors);
secondResors4b50_CR24=secondResors([1:1:mincolsize],[1:1:size(secondResors,2)]);
firstResors4b50_CR24=firstResors([1:1:mincolsize],[1:1:size(firstResors,2)]);
gatherInfo4b50_CR24=gatherInfo([1:1:mincolsize],[1:1:size(gatherInfo,2)]);
finalHead50_CR24=finalHead;

[ totalNum totalInfo MapCR FidRatio  QuantiBits Density BestRound gatherInfo supNum firstEnergy SecondEnergy ... 
    firstResors secondResors headx finalHead] = parseResorsToMatrixByHead_4b('data/ULSA4b4_All_N75_BW180.0PW0.001_FR0.95_r500.0.txt');
mincolsize=returnColNonZeroSize(firstResors);
secondResors4b75_CR24=secondResors([1:1:mincolsize],[1:1:size(secondResors,2)]);
firstResors4b75_CR24=firstResors([1:1:mincolsize],[1:1:size(firstResors,2)]);
gatherInfo4b75_CR24=gatherInfo([1:1:mincolsize],[1:1:size(gatherInfo,2)]);
finalHead75_CR24=finalHead;

[ totalNum totalInfo MapCR FidRatio  QuantiBits Density BestRound gatherInfo supNum firstEnergy SecondEnergy ... 
    firstResors secondResors headx finalHead] = parseResorsToMatrixByHead_4b('data/ULSA4b4_All_N100_BW180.0PW0.001_FR0.95_r500.0.txt');
mincolsize=returnColNonZeroSize(firstResors);
secondResors4b100_CR24=secondResors([1:1:mincolsize],[1:1:size(secondResors,2)]);
firstResors4b100_CR24=firstResors([1:1:mincolsize],[1:1:size(firstResors,2)]);
gatherInfo4b100_CR24=gatherInfo([1:1:mincolsize],[1:1:size(gatherInfo,2)]);
finalHead100_CR24=finalHead;

[ totalNum totalInfo MapCR FidRatio  QuantiBits Density BestRound gatherInfo supNum firstEnergy SecondEnergy ... 
    firstResors secondResors headx finalHead] = parseResorsToMatrixByHead_4b('data/ULSA4b4_All_N125_BW180.0PW0.001_FR0.95_r500.0.txt');
mincolsize=returnColNonZeroSize(firstResors);
secondResors4b125_CR24=secondResors([1:1:mincolsize],[1:1:size(secondResors,2)]);
firstResors4b125_CR24=firstResors([1:1:mincolsize],[1:1:size(firstResors,2)]);
gatherInfo4b125_CR24=gatherInfo([1:1:mincolsize],[1:1:size(gatherInfo,2)]);
finalHead125_CR24=finalHead;

[ totalNum totalInfo MapCR FidRatio  QuantiBits Density BestRound gatherInfo supNum firstEnergy SecondEnergy ... 
    firstResors secondResors headx finalHead] = parseResorsToMatrixByHead_4b('data/ULSA4b4_All_N150_BW180.0PW0.001_FR0.95_r500.0.txt');
mincolsize=returnColNonZeroSize(firstResors);
secondResors4b150_CR24=secondResors([1:1:mincolsize],[1:1:size(secondResors,2)]);
firstResors4b150_CR24=firstResors([1:1:mincolsize],[1:1:size(firstResors,2)]);
gatherInfo4b150_CR24=gatherInfo([1:1:mincolsize],[1:1:size(gatherInfo,2)]);
finalHead150_CR24=finalHead;


[ totalNum totalInfo MapCR FidRatio  QuantiBits Density BestRound gatherInfo supNum firstEnergy SecondEnergy ... 
    firstResors secondResors headx finalHead] = parseResorsToMatrixByHead_4b('data/ULSA4b4_All_N195_BW180.0PW0.001_FR0.95_r500.0.txt');
mincolsize=returnColNonZeroSize(firstResors);
secondResors4b195_CR24=secondResors([1:1:mincolsize],[1:1:size(secondResors,2)]);
firstResors4b195_CR24=firstResors([1:1:mincolsize],[1:1:size(firstResors,2)]);
gatherInfo4b195_CR24=gatherInfo([1:1:mincolsize],[1:1:size(gatherInfo,2)]);
finalHead195_CR24=finalHead;


[ totalNum totalInfo MapCR FidRatio  QuantiBits Density BestRound gatherInfo supNum firstEnergy SecondEnergy ... 
    firstResors secondResors headx finalHead] = parseResorsToMatrixByHead_4b('data/ULSA4b4_All_N50_BW180.0PW0.001_FR0.95_r500.36.txt');
mincolsize=returnColNonZeroSize(firstResors);
secondResors4b50_CR36=secondResors([1:1:mincolsize],[1:1:size(secondResors,2)]);
firstResors4b50_CR36=firstResors([1:1:mincolsize],[1:1:size(firstResors,2)]);
gatherInfo4b50_CR36=gatherInfo([1:1:mincolsize],[1:1:size(gatherInfo,2)]);
finalHead50_CR36=finalHead;

[ totalNum totalInfo MapCR FidRatio  QuantiBits Density BestRound gatherInfo supNum firstEnergy SecondEnergy ... 
    firstResors secondResors headx finalHead] = parseResorsToMatrixByHead_4b('data/ULSA4b4_All_N75_BW180.0PW0.001_FR0.95_r500.36.txt');
mincolsize=returnColNonZeroSize(firstResors);
secondResors4b75_CR36=secondResors([1:1:mincolsize],[1:1:size(secondResors,2)]);
firstResors4b75_CR36=firstResors([1:1:mincolsize],[1:1:size(firstResors,2)]);
gatherInfo4b75_CR36=gatherInfo([1:1:mincolsize],[1:1:size(gatherInfo,2)]);
finalHead75_CR36=finalHead;

[ totalNum totalInfo MapCR FidRatio  QuantiBits Density BestRound gatherInfo supNum firstEnergy SecondEnergy ... 
    firstResors secondResors headx finalHead] = parseResorsToMatrixByHead_4b('data/ULSA4b4_All_N100_BW180.0PW0.001_FR0.95_r500.36.txt');
mincolsize=returnColNonZeroSize(firstResors);
secondResors4b100_CR36=secondResors([1:1:mincolsize],[1:1:size(secondResors,2)]);
firstResors4b100_CR36=firstResors([1:1:mincolsize],[1:1:size(firstResors,2)]);
gatherInfo4b100_CR36=gatherInfo([1:1:mincolsize],[1:1:size(gatherInfo,2)]);
finalHead100_CR36=finalHead;

[ totalNum totalInfo MapCR FidRatio  QuantiBits Density BestRound gatherInfo supNum firstEnergy SecondEnergy ... 
    firstResors secondResors headx finalHead] = parseResorsToMatrixByHead_4b('data/ULSA4b4_All_N125_BW180.0PW0.001_FR0.95_r500.36.txt');
mincolsize=returnColNonZeroSize(firstResors);
secondResors4b125_CR36=secondResors([1:1:mincolsize],[1:1:size(secondResors,2)]);
firstResors4b125_CR36=firstResors([1:1:mincolsize],[1:1:size(firstResors,2)]);
gatherInfo4b125_CR36=gatherInfo([1:1:mincolsize],[1:1:size(gatherInfo,2)]);
finalHead125_CR36=finalHead;

[ totalNum totalInfo MapCR FidRatio  QuantiBits Density BestRound gatherInfo supNum firstEnergy SecondEnergy ... 
    firstResors secondResors headx finalHead] = parseResorsToMatrixByHead_4b('data/ULSA4b4_All_N150_BW180.0PW0.001_FR0.95_r500.36.txt');
mincolsize=returnColNonZeroSize(firstResors);
secondResors4b150_CR36=secondResors([1:1:mincolsize],[1:1:size(secondResors,2)]);
firstResors4b150_CR36=firstResors([1:1:mincolsize],[1:1:size(firstResors,2)]);
gatherInfo4b150_CR36=gatherInfo([1:1:mincolsize],[1:1:size(gatherInfo,2)]);
finalHead150_CR36=finalHead;


[ totalNum totalInfo MapCR FidRatio  QuantiBits Density BestRound gatherInfo supNum firstEnergy SecondEnergy ... 
    firstResors secondResors headx finalHead] = parseResorsToMatrixByHead_4b('data/ULSA4b4_All_N195_BW180.0PW0.001_FR0.95_r500.36.txt');
mincolsize=returnColNonZeroSize(firstResors);
secondResors4b195_CR36=secondResors([1:1:mincolsize],[1:1:size(secondResors,2)]);
firstResors4b195_CR36=firstResors([1:1:mincolsize],[1:1:size(firstResors,2)]);
gatherInfo4b195_CR36=gatherInfo([1:1:mincolsize],[1:1:size(gatherInfo,2)]);
finalHead195_CR36=finalHead;


%plot(headx,mean(finalHead195),'*-b','LineWidth',1.5,'DisplayName','Head\leqm (195 Machines)','MarkerSize',10);hold on;
%plot(headx,mean(finalHead150),'^--r','LineWidth',1.5,'DisplayName','Head\leqm (150 Machines)','MarkerSize',10);hold on;
%plot(headx,mean(finalHead125),'^--r','LineWidth',1.5,'DisplayName','Head\leqm (125 Machines)','MarkerSize',10);hold on;
%plot(headx,mean(finalHead100),'^--r','LineWidth',1.5,'DisplayName','Head\leqm (100 Machines)','MarkerSize',10);hold on;
%plot(headx,mean(finalHead75),'kx-.','LineWidth',1.5,'DisplayName','Head\leqm (75 Machines)','MarkerSize',10);hold on;
%plot(headx,mean(finalHead50),'kx-.','LineWidth',1.5,'DisplayName','Head\leqm (50 Machines)','MarkerSize',10);hold on;

machineNum = [50 75 100 125 150 195];
idx = 5;
headNum = [ mean(finalHead50(:,idx)) mean(finalHead75(:,idx)) mean(finalHead100(:,idx)) mean(finalHead125(:,idx)) mean(finalHead150(:,idx)) mean(finalHead195(:,idx))];
idx24 = 5;
machineNum24 = [50 75 100 125 150 195];
headNum24 = [ mean(finalHead50_CR24(:,idx24)) mean(finalHead75_CR24(:,idx24)) mean(finalHead100_CR24(:,idx24)) mean(finalHead125_CR24(:,idx24)) mean(finalHead150_CR24(:,idx24))  mean(finalHead195_CR24(:,idx24))];
idx36 = 5;
machineNum36 = [50 75 100 125 150 195];
headNum36 = [ mean(finalHead50_CR36(:,idx36)) mean(finalHead75_CR36(:,idx36)) mean(finalHead100_CR36(:,idx36)) mean(finalHead125_CR36(:,idx36)) mean(finalHead150_CR36(:,idx36))  mean(finalHead195_CR36(:,idx36))];
plot(machineNum, headNum,'*-b','LineWidth',1.5,'DisplayName','Increase of HeadNumber CR = 0.48','MarkerSize',10); hold on;
%plot(machineNum, headNum36,'kx-.','LineWidth',1.5,'DisplayName','Increase of HeadNumber CR = 0.36','MarkerSize',10); hold on;
plot(machineNum, headNum24,'*-r','LineWidth',1.5,'DisplayName','Increase of HeadNumber CR = 0.24','MarkerSize',10); hold on;


%plot(headx,mean(finalHead),'^-.','LineWidth',1.5,'DisplayName','Head\leqm','MarkerSize',10);hold on;

%plot (headx,headx,'kx-','LineWidth',1.5,'DisplayName','Head=m');

title('180 Khz;\lambda=0.95;\eta=0.48');
ylabel('Average Support Head Number');
xlabel({'Number of Machine','Two-tier Data Gathering'});
legend('show');
grid on;
%axis([3 23 4 23 ]);
