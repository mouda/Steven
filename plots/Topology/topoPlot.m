clear; clear all;

dmp = ...
    dlmread('../../runSimulation/runHN/mapFile/mapFile_uni_50_r500/mapFile_uniR500_N50_1.txt');

numMachine = dmp(1,1);
radius = dmp(1,2);
summarizeOffset = 1;
if size(dmp,1)-summarizeOffset ~= numMachine
    fprintf('wrang number');
    exit;
end
x = dmp(1+summarizeOffset : numMachine + summarizeOffset , 1);
y = dmp(1+summarizeOffset : numMachine + summarizeOffset , 2);

circle(0,0,radius); hold on;
scatter(x,y); 
grid on;