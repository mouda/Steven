function [ InfoRatio NumRatio totalResors firstResors SecResors Round ] = parse_EachIteration( input_args )
%PARSE_EACHITERATION Summary of this function goes here
%   Detailed explanation goes here
A=dlmread(input_args);
InfoRatio=A(:,1);
NumRatio=A(:,2);
totalResors=A(:,3);
firstResors=A(:,4);
SecResors=A(:,5);
Round=A(:,6);
end

