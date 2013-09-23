%Create By: Steven
%Thanks: Alvin and Jack @Tonic
%Purpose:
%   draw the Topology and create interference contour

function []=showSINRContour(mapPath,StructurePath)

%Marco:---------
rawToDB = @(raw) 10*log10(raw);
dBToRaw = @(dB) 10.^(dB/10);
CombineDB = @(dB1,dB2) 10*log10(10.^(dB1/10) + 10.^(dB2/10));
%---------------
noise = 180 * 10^3 *10^-19;
[ totalNodes radius x y Gij Gib ] = parse_Map_v1( mapPath );
[ maxChNum powerMax C2WT payoffs SAFac TimeMs1st TimeMs2nd EnergyJ1st EnergyJ2nd clusterStru headList TxPower  ] = parse_Stru_v1 ...
    (StructurePath,totalNodes);

figure;
circle(0,0,radius);hold on;
hold on;
% Create plot
plot(0,0,'MarkerFaceColor',[0 0 0],'MarkerSize',10,'Marker','^','LineStyle','none', 'DisplayName','Base Station');
axis equal;
% Create xlabel
xlabel({'x-axis (m)', 'For Sensor Field'});
% Create ylabel
ylabel({'y-axis (m)'});

densityScale = 1000;%km square
density = totalNodes / (pi*(radius/densityScale)^2) ;


B = find(sum(clusterStru));
supSet = setdiff(B, headList);

numOfTx = max(size(headList)); % plus base station
plot(x(supSet),y(supSet),'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1],'Marker','o',...
   'LineStyle','none', 'MarkerSize', 7, 'DisplayName','Supported machine');
plot(x(headList),y(headList),'MarkerFaceColor',[0.3 0.3 0.3],'MarkerEdgeColor',[0.3 0.3 0.3],'Marker','s',...
   'LineStyle','none','MarkerSize',5, 'DisplayName','Cluster Head');

plot(x,y,'MarkerFaceColor',[0.3 0.3 0.3],'MarkerEdgeColor',[0.3 0.3 0.3],...
    'MarkerSize',7,...
    'Marker','x',...
    'LineStyle','none', 'DisplayName','Machine');



TxPower = rawToDB(TxPower * 1000); % in dBm

BS_loc = [0 0];

step = 1.0;
[X,Y] = meshgrid(-radius:step:radius,-radius:step:radius);

resol = size(-radius:step:radius, 2);

%InterAll: a 3-D cube, each element means the most strongest interference 
%   from each cluster j. e.g.(1,1,j) the strongest interference from
%   cluster j in grid point 1,1;
InterAll=-inf.*ones(resol,resol,maxChNum); 
for i=1:maxChNum
    oneCluX=x(find(clusterStru(i,:)));
    oneCluY=y(find(clusterStru(i,:)));
    oneCluPower = TxPower(find(clusterStru(i,:)));
        for j=1:size(oneCluX) 
        tmpX_MA = oneCluX(j).*ones(resol, resol);
        tmpY_MA = oneCluY(j).*ones(resol, resol);
        tmpP_MA = oneCluPower(j).*ones(resol, resol);
        RxPower1 = prop_MA(tmpP_MA,tmpX_MA,tmpY_MA,X,Y);
        InterAll(:,:,i)=max(InterAll(:,:,i),RxPower1);
        end
end
%Change from db to power
InterAll=dBToRaw(InterAll);

for i = 1:maxChNum
    interf=sum(InterAll,3) - InterAll(:,:,i);
    oneCluX=x(find(clusterStru(i,:)));
    oneCluY=y(find(clusterStru(i,:)));
    oneCluPower = TxPower(find(clusterStru(i,:)));
    sizeOfClu = size(oneCluX,1)-1;
    
    for j=1:size(oneCluX)    
        tmpX_MA = oneCluX(j).*ones(resol, resol);
        tmpY_MA = oneCluY(j).*ones(resol, resol);
        tmpP_MA = oneCluPower(j).*ones(resol, resol);
        RxPower1 = prop_MA(tmpP_MA,tmpX_MA,tmpY_MA,X,Y);
        RxPower1 = dBToRaw(RxPower1);
        SINR = RxPower1./(interf+noise);
        Delta=0.05;
        rate = log2(1+SINR);
        %figure;
        contour(X,Y,rate,sizeOfClu*C2WT);hold on;
    end
end

end