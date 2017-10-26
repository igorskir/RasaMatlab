%% Load the distributions 
isNormalize = 0;
if isNormalize == 0
    sample = struct2mat(dataTraining);
else
    sample = dataTrainingNorm;
end
featureArr = sample(:,1);
featureArr = sample(:,10);
vars.distributionLoading = {'isNormalize', 'sample'};
clear(vars.distributionLoading{:});
%%
[numBins,edges] = histcounts(featureArr, 'BinMethod', 'fd');
numBinsCatheter = zeros(1, size(numBins,2));
numBinsTissue = zeros(1, size(numBins,2));
numRows = size(featureArr, 1);
numCols = size(featureArr, 2);
intervalStep = edges(2) - edges(1);
% intervalStep = round(edges(end)/size(numBins,2));
for i = 1:numRows
    bin = floor((featureArr(i,2) - edges(1))/intervalStep) + 1; 
        if featureArr(i,1) == 1 
            numBinsCatheter(1, bin) = numBinsCatheter(1, bin) + 1;
        else
            numBinsTissue(1, bin) = numBinsTissue(1, bin) + 1;
        end  
end
vars.separation = {'numRows', 'numCols', 'i', 'bin', 'numBins'};
clear(vars.separation{:});

%% Normalization
numBinsCatheterNorm = numBinsCatheter/size(dataCatheter,1);
numBinsTissueNorm = numBinsTissue/size(dataTissue,1);
catheterSum = sum(numBinsCatheterNorm);
catheterSum = sum(numBinsTissueNorm);
vars.normalization = {'catheterSum', 'catheterSum'};
clear(vars.normalization{:});

%% Get the score
overallDelta = 0;
for i = 1:size(numBinsCatheterNorm,2)
    currentDelta = numBinsCatheterNorm(1,i) - numBinsTissueNorm(1,i);
    if numBinsCatheterNorm(1,i) ~= 0 && currentDelta > 0
        overallDelta = overallDelta + currentDelta; 
    end    
end

plot(numBinsCatheterNorm)
hold on
plot(numBinsTissueNorm)

vars.scoreObtaining = {'i', 'currentDelta', 'numBinsCatheterNorm', 'numBinsTissueNorm'};
clear(vars.scoreObtaining{:});
