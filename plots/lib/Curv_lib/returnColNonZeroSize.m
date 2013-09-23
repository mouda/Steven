function [ output_args ] = returnColNonZeroSize( input_args )
%output_args = minimum number of non zero element of allcolumns
%input_args has to be a matrix
output_args=100000;%Arbitarary large number
for i=1:size(input_args,2)
    if(length(find(input_args(:,i)))<output_args)
        output_args=length(find(input_args(:,i)));
    end
end
end

