function [ maxChNum powerMax C2W payoffs SAFac TimeMs1st TimeMs2nd EnergyJ1st EnergyJ2nd clusterStru headList TxPower  ] = parse_Stru_v1( input_paths,expTotalNodes )
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
  C2W = fscanf(fid,'%e',1);
  payoffs=fscanf(fid,'%d',1);
  SAFac=fscanf(fid,'%d',1);  
  TimeMs1st= fscanf(fid,'%e',1);
  TimeMs2nd = fscanf(fid,'%e',1);%the record is wrong
  EnergyJ1st = fscanf(fid,'%e',1);
  EnergyJ2nd  = fscanf(fid,'%e',1);  
else
  error('Structure file does not match topology');
end

clusterStru = zeros(maxChNum,expTotalNodes);
for i=1:maxChNum
  clusterStru(i,:) = fscanf(fid,'%d',expTotalNodes);
end
headList = fscanf(fid,'%d',maxChNum);
TxPower = fscanf(fid,'%f',expTotalNodes);
fclose(fid);
end

