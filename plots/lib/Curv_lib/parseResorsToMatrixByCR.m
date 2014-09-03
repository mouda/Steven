function [ totalNum totalInfo MapCR FidRatio  QuantiBits Density BestRound gatherInfo supNum firstEnergy SecondEnergy firstResors secondResors CRx] = parseResorsToMatrixByCR( input_args )
data=dlmread(input_args);
CR = data(:,4);
[CV CI]=sort(CR);
CRx=unique(CR);   
for i =1:length(CRx)
    target=CRx(i);
    arr=find(data(:,4)==target);
    for j=1:length(arr)
        totalNum(j,i)= data(arr(j),1);
        totalInfo(j,i)= data(arr(j),3);
        MapCR(j,i) = data(arr(j),4);
        FidRatio(j,i) = data(arr(j),5);
        QuantiBits(j,i) = data(arr(j),6);
        Density(j,i) = data(arr(j),7);
        BestRound(j,i) = data(arr(j),8);
        gatherInfo(j,i)= data(arr(j),9);
        supNum(j,i)= data(arr(j),10);
        firstEnergy(j,i)=data(arr(j),11);
        SecondEnergy(j,i)=data(arr(j),12);
        firstResors(j,i)=data(arr(j),13);
        secondResors(j,i)=data(arr(j),14);
    end
end
end

