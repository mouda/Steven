clear all;
radius = 150;%(meter)
area = pi * radius^2 ;
% use this to generate the mapfile
expA=[10 20 30 40 50];  

for i =1:length(expA)
    display start
    exp = expA(i);
    for runTime = 1:50
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
        folderName = sprintf('mapFile_uni_%d_r%d',exp, radius);
        if(~exist(folderName, 'file'))
            mkdir(folderName);
        end
        fid = fopen(sprintf('%s/mapFile_uniR%d_N%d_%d.txt', folderName,radius, exp, runTime), 'wt');
        fprintf(fid, '%d %d\n', exp, radius);
        fprintf(fid,'%f\t%f\n',[x; y]);
        fclose(fid);
    end
    display done;
    close all;
end
