data = dlmread('data/estimation');
data11 = dlmread('data/estimate11');
data19 = dlmread('data/estimate19');
dataLength = size(data,1);
%x = zeros(dataLength,1);
%y = zeros(dataLength,1);
x = [];
y = [];
%for i = 1:size(data,1)
%    if data(i,3) - data(i,1) < 0.0
%         x = [x (data(i,2) - data(i,1))/data(i,1)];
%    y = [y (data(i,3) - data(i,1))/data(i,1)];
%    end
   
%end

mean(data(:,4))
%scatter(x,y,'*');
%axis([-0.02 0.08 -0.05 1 ])
%grid on;
%lsline;

target = [ mean(mean(data11(:,4))) mean(mean(data19(:,4))) mean(mean(data(:,4)))];
index=[1 2 3 ];
Obar = bar(index, target);
grid on;
%set(Obar(1),'FaceColor',[1 0.35 0.35]);
%set(Obar(2),'FaceColor',[0 1 0.7]);
%set(Obar(1),'FaceColor',[1 1 0.9]);
ylabel('Error percentage (%)');

set(gca,'XTick',index);
xlabel({'Maximum Cluster Head'});  
set(gca,'XTickLabel',{'11';'19';'29'});
