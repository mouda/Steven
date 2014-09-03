clear all;
radius = 100;%(meter)
area = pi * radius^2 ;

expA=[100 150 200 250 300];  
for i =1:5
    exp = expA(i);
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
        folderName = sprintf('map_uni_den');
        if(~exist(folderName, 'file'))
            mkdir(folderName);
        end
        fid = fopen(sprintf('%s/mapFile_uniR%d_N%d_1.txt', folderName,radius, exp), 'wt');
        fprintf(fid, '%d\t%d\n', exp,radius);
        fprintf(fid,'%f\t%f\n',[x; y]);
        fclose(fid);
    end

    close all;
end
