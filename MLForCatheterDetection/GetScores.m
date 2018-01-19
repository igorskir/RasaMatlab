%% Estimate scores of PDF distribution (SBFS)
clc; clear;
addpath('MAT files');
load('data (not separ.).mat');
binCounts = 10:2:70;
featureCount = size(dataCatheter,2);
scoresSBFS = zeros(numel(binCounts),featureCount);
for i = 1:numel(binCounts)
    scoresSBFS(i,1) = binCounts(i);
    for j = 2:featureCount
        isDiscrete = IsDiscrete(dataCatheter(:,j));
        scoresSBFS(i,j) = ShowPDFDifference(dataTissue(:,j), dataCatheter(:,j), binCounts(i), binCounts(i), isDiscrete);
%         scores(i,j) = ShowPDFDifference(dataTissue(:,j), dataCatheter(:,j), binCounts(i), binCounts(i), isDiscrete, 'sPlot');
%         close
    end
end

%% Estimate scores of overlapping (OFS)
epsilonCounts = 0.01:0.01:0.10;
scoresOFS = zeros(1, featureCount);
for i = 1:numel(epsilonCounts)
    scoresOFS(i,1) = epsilonCounts(i);
    for j = 2:featureCount
        scoresOFS(i,j) = GetOverlapRate(dataTissue(:,j), dataCatheter(:,j), epsilonCounts(i));
    end
end
