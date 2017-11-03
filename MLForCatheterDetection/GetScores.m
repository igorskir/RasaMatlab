clc;
clear;
addpath('MAT files');
load data;
binCounts=10:10:100;
featureCount=size(dataCatheter,2);
scores=zeros(numel(binCounts),featureCount);
for i=1:numel(binCounts)
    scores(i,1)=binCounts(i);
    for j=2:featureCount
        isDiscrete=IsDiscrete(dataCatheter(:,j));
        scores(i,j)=ShowPDFDifference(dataTissue(:,j),dataCatheter(:,j),binCounts(i),binCounts(i),isDiscrete);
    end
end
