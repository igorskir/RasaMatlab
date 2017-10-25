%% Load the distributions 

isNormalize = 0;
if isNormalize == 0
    sample = struct2mat(dataTraining);
else
    sample = dataTrainingNorm;
end
[numBins,edges] = histcounts(sample, 'BinMethod', 'auto');


%%
numBinsCatheter = zeros(1, size(numBins,2));
numBinsTissue = zeros(1, size(numBins,2));
numRows = size(sample, 1);
numCols = size(sample, 2);
for i = 1:numRows
    novayaPeremennaya = floor((sample(i,2) - edges(1))/(edges(2) - edges(1))) + 1; 
        if sample(i,1) == 1 
            numBinsCatheter(1, novayaPeremennaya) = numBinsCatheter(1, novayaPeremennaya) + 1;
        else
            numBinsTissue(1, novayaPeremennaya) = numBinsTissue(1, novayaPeremennaya) + 1;
        end  
end
numCathCases = sum(numBinsCatheter);
numTissueCases = sum(numBinsTissue);

numBinsCatheter = numBinsCatheter/(numCathCases*(edges(2) - edges(1)));
numBinsTissue = numBinsTissue/(numTissueCases*(edges(2) - edges(1)));

plot(numBinsTissue);
hold on
plot(numBinsCatheter);

numBinsCatheter(numBinsCatheter == 0) = [];
numBinsTissue(numBinsTissue == 0) = [];
z = numBinsCatheter./(numBinsCatheter + numBinsTissue);

featureScore = mean(z); 






