%V2 is for parsing structure that second tier resource is a variable;
function [ maxChNum powerMax C2W payoffs SAFac TimeMs1st TimeMs2nd EnergyJ1st EnergyJ2nd clusterStru headList TxPower  ] = parse_Stru_v2( input_paths,expTotalNodes )
%PARSE_STRU_V1 Summary of this function goes here
%  parse the structure before and after 2013/05/14
[fid message] = fopen(input_paths,'r');
if fid ==-1
  disp(message);
  error('Read Structure Fail');
end
if expTotalNodes == fscanf(fid,'%d',1);
  maxChNum = fscanf (fid,'%d',1);
  powerMax = fscanf(fid,'%e',1);
  indEntropy = fscanf(fid,'%e',1);
  payoffs=fscanf(fid,'%e',1);
  SAFac=fscanf(fid,'%d',1);  
  TimeMs1st= fscanf(fid,'%e',1);
  TimeMs2nd = fscanf(fid,'%e',1);
  EnergyJ1st = fscanf(fid,'%e',1);
  EnergyJ2nd  = fscanf(fid,'%e',1);
else
  error('Structure file does not match topology');
end
C2W=indEntropy/TimeMs2nd/180;%hardcode 180k HZ
clusterStru = zeros(maxChNum,expTotalNodes);
for i=1:maxChNum
  clusterStru(i,:) = fscanf(fid,'%d',expTotalNodes);
end
headList = fscanf(fid,'%d',maxChNum);
TxPower = fscanf(fid,'%f',expTotalNodes);
fclose(fid);
end

