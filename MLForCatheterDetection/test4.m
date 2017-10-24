clear all;
load('a3.mat')
a3Target = a3(1:end, 1);
a3(:,1) = [];
a3Input = a3;
%%
a3InputNorm = robustNormalization(a3Input, 'quantile', 1);
%%
[numRows, numCols] = size(a3InputNorm);
a3InputNorm = a3InputNorm';
netSimResults = zeros(numRows,1);
tic;
for i = 1:numRows
    netSimResults(i,:) = round(net(a3InputNorm(:,i)));
end
load('a3.mat')