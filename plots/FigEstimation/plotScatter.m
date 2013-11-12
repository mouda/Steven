data = dlmread('data/estimation');

dataLength = size(data,1);
x = zeros(dataLength,1);
y = zeros(dataLength,1);
for i = 1:size(data,1)
    x(i) = (data(i,2) - data(i,1))/data(i,1);
    y(i) = (data(i,3) - data(i,1))/data(i,1);
end

scatter(x,y,'*');
axis([-0.02 0.08 -0.05 1 ])
grid on;