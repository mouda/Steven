function createfigure(XMatrix1, YMatrix1, EMatrix1)
%CREATEFIGURE(XMATRIX1,YMATRIX1,EMATRIX1)
%  XMATRIX1:  errorbar x matrix
%  YMATRIX1:  errorbar y matrix
%  EMATRIX1:  errorbar e matrix

%  Auto-generated by MATLAB on 21-Nov-2013 13:26:13

% Create figure
figure1 = figure('XVisual',...
    '0xaf (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)');

% Create axes
axes1 = axes('Parent',figure1);
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');

% Create multiple error bars using matrix input to errorbar
errorbar1 = errorbar(XMatrix1,YMatrix1,EMatrix1,'MarkerSize',10,...
    'LineWidth',1.5);
set(errorbar1(1),'Marker','^','DisplayName','Two-Tier DC SA',...
    'Color',[0 0 1]);
set(errorbar1(2),'Marker','*','DisplayName','Distance kmeans (m=11)',...
    'Color',[0.5 0.5 1]);
set(errorbar1(3),'Marker','^','LineStyle','--',...
    'DisplayName','DC Kmeans (m=11)',...
    'Color',[0 0 0]);
set(errorbar1(4),'Marker','o','LineStyle','-.',...
    'DisplayName','MC Kmeans (m=11)',...
    'Color',[1 0 0]);

% Create xlabel
xlabel({'Compression Ratio (\eta)','Two-Tier Data Gathering'});

% Create ylabel
ylabel('Resource Usage Ratio (R)');

% Create title
title('195 Machines;180 Khz;\lambda=0.95');

% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.533333333333333 0.748611111111111 0.37114583333334 0.175277777777778]);

