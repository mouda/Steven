function [Rate SINR Inter ] = calSINR_4All( maxChNum, totalNodes, clusterStru,  clusterSize,  headName, TxPower,  Gij  )
%CALSINR_4ALL Summary of this function goes here
%   Detailed explanation goes here
BW = 180*1e3;
n0 = BW*1e-19;
SINR = zeros(1,totalNodes);
for i = 1:maxChNum
   %find the Interference from each node first;
   
  if(headName(i)==0)continue;end
   accu = n0;
   for j=1:maxChNum
      if (i~=j)
         maxInter = -1;
        for k=1:totalNodes
            if(clusterStru(j,k)==1)
                tempInter = TxPower(k)*Gij(headName(i),k);
                if(tempInter>maxInter)
                   maxInter = tempInter;
                end
            end
        end
        accu = accu + maxInter;
      end
   end
   Inter(i) = accu;
   %after get interference we compute SINR for each node
   for j=1:totalNodes
     if(clusterStru(i,j)==1&&headName(i)==j)
         SINR(j)=0;
     elseif(clusterStru(i,j)==1)
         SINR(j) = TxPower(j)*Gij(headName(i),j)/accu;
         Rate(j) = BW * log2(1+TxPower(j)*Gij(headName(i),j)/accu)/(clusterSize(i)-1);
     end
   end
   
end

end

