function [result] = ParseToPart(fileName, k, nrep)
weightMatrix = dlmread(fileName);
result = grPartition(weightMatrix, k,nrep);
str = sprintf('%s_part_%d',fileName,k);
result = result -1;
dlmwrite(str, result);
end


