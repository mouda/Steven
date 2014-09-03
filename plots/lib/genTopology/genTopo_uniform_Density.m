clear all;
aryDensity = [2840 3500 4000 4500 5000 5500 6000 6500 6694 7500 8500 9000 9500 10000];
radius = 500;%(meter)
area = pi * radius^2 ;

for i =1:length(aryDensity)
    aryTotalNodes(i) = floor( area *(1e-6) *aryDensity(i));
end
  
for i =1:length(aryDensity)
    exp = aryTotalNodes(i);
    for runTime = 1:1
        x = zeros(1, exp);
        y = zeros(1, exp);
        for j=1:exp
            x(j) = radius;
            y(j) = radius;
            while (x(j)^2+y(j)^2 > radius^2) 
                x(j) = rand() * 2*radius - radius;
                y(j) = rand() * 2*radius - radius;
            end
        end
       
       %{ 
        figure;
        hold on;
        grid on;
        plot(x,y,'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7],...
                  'MarkerSize',4,...
                  'Marker','o',...
                  'LineStyle','none');
        axis equal;
        circle(radius, 0, 0, 'k', 10000);
        xlabel({'x-axis (m)'});
        ylabel({'y-axis (m)'});
        xlim([-100 100]);
        ylim([-100 100]);
        set(gca,'XTick', -radius:20:radius,...
                'YTick', -radius:20:radius);
        %}
        folderName = sprintf('map_uniform_Dsensity', exp,i);
        if(~exist(folderName, 'file'))
            mkdir(folderName);
        end
        fid = fopen(sprintf('%s/mapFile_uniform_N%d_Den%d_%d.txt', folderName, exp, i, runTime), 'wt');
        fprintf(fid, '%d\n', exp);
        fprintf(fid,'%f\t%f\n',[x; y]);
        fclose(fid);
    end

    close all;
end
%}