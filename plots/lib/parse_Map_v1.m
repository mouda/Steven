function [ totalNodes radius x y Gij Gib ] = parse_Map_v1( input_path )
%PARSE_MAP_V1 Summary of this function goes here
%   Detailed explanation goes here
%fid = fopen('mapFile/mapFile_uni_195_r500/mapFile_uniR500_N195_2.txt', 'r');

rawToDB = @(raw) 10*log10(raw);
dBToRaw = @(dB) 10.^(dB/10);
CombineDB = @(dB1,dB2) 10*log10(10.^(dB1/10) + 10.^(dB2/10));
[fid message] = fopen(input_path, 'r');
if fid ==-1
  disp(message);
end
C = textscan(fid, '%f %f', 1);
totalNodes = C{1};
radius=C{2};
C = textscan(fid, '%f %f', totalNodes);
x = C{1};
y = C{2};
[X1, X2] = meshgrid(x(1:totalNodes));
[Y1, Y2] = meshgrid(y(1:totalNodes));
X_sq = (X1-X2).^2;
Y_sq = (Y1-Y2).^2;
D2 = X_sq + Y_sq;
Gij = dBToRaw(-(131.1 + 42.81*log10( (D2.^0.5)/1000 )));
Dbs=x.^2+y.^2;
Gib = dBToRaw(-(131.1 + 42.81*log10( (Dbs.^0.5)/1000 )));
end

