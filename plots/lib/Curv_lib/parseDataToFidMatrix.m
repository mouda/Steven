function [ FidMatrix Resorsx headx] = parseDataToFidMatrix( input_args )

data=dlmread(input_args);
wholeSysInfo=data(1,3);
head = data(:,2);
Resors = data(:,14);
[HV HI]=sort(head);
headx=unique(HV);
[RV RI]=sort(Resors);
Resorsx=unique(Resors);    

for i =1:length(headx)
    target=headx(i);
    arr=find(data(:,2)==target);
    for j=1:length(Resorsx)
        Rtarget=Resorsx(j);
        Rarr=find(data(arr,14)==Rtarget); 
        FidMatrix(i,j)=mean(data(arr(Rarr),9)./wholeSysInfo);
    end
end
end

