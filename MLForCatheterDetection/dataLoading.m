%% Initial state
clear all; close all; clc;
set(0, 'DefaultFigureWindowStyle', 'normal');
currentFolder = pwd;
addpath(genpath(pwd));
isGetDistance = 0;
isNormalize = 0;
%% Load the data            
dataTraining = struct('Presence', [], ...
                      'Area', [], ...
                      'ConvexArea', [], ...
                      'Perimeter', [],  ...
                      'Eccentricity', [], ...
                      'Solidity', [], ...
                      'Extent', [], ...
                      'EquivDiameter', [], ...
                      'MaxIntensity', [], ...
                      'MeanIntensity', [], ...
                      'MinIntensity', [], ...
                      'Variance', [], ...
                      'StandardDeviation' , [], ...
                      'Contrast', [], ...
                      'Correlation', [], ...
                      'Energy', [], ...
                      'Homogeneity', [], ...
                      'MajorAxisLength', [], ... 
                      'MinorAxisLength', []);

% Load training data
% trainingFolder = strcat(currentFolder, '\Labeled data\LV Catheter 07 v.3 (18 features)\Training');
trainingFolder = strcat(currentFolder, '\Labeled data\LV Catheter 07 v.4 (18 features)\Training');
cd(trainingFolder);
s = what;
matfiles = s.mat;
for i = 1:numel(matfiles)
    tempTrain = load(char(matfiles(i)));
    tempTrain = tempTrain.featuresAll;
    dataTraining = [tempTrain, dataTraining];
end
dataTraining(end) = [];
cd ..\..\..\
vars.loadTrainData = {'currentFolder', 'trainingFolder' , 'i', 's', ...
                      'matfiles', 'tempTrain'};
clear(vars.loadTrainData{:});

%% Removing NaNs
fn = fieldnames(dataTraining);
for i = 1:numel(fn) 
    dataTraining = dataTraining(~isnan([dataTraining.(fn{i})]));
end
vars.removingNaNs = {'fn', 'i'};
clear(vars.removingNaNs{:});

%% Get normalized data
if isNormalize == 1
    dataTrainingNorm = robustNormalization(dataTraining, 'quantile'); % quartile, std, mean, linear
    numFields = numel(fieldnames(dataTraining));
    dataTissueNorm = zeros(1, numFields, 'double');
    dataCatheterNorm = zeros(1, numFields, 'double');
    for i = 1:size(dataTrainingNorm, 1)
       if dataTrainingNorm(i,1) == 1
            dataCatheterNorm(end+1,:) = dataTrainingNorm(i,:); 
       else 
            dataTissueNorm(end+1,:) = dataTrainingNorm(i,:);  
       end   
    end
    netTrainTargetsNorm = dataTrainingNorm(:,1);
    netTrainInputsNorm = dataTrainingNorm;
    netTrainInputsNorm(:,1) = [];
    dataCatheterNorm(1,:) = [];
    dataTissueNorm(1,:) = [];
    vars.normalization = {'numFields', 'dataTemp', 'i'};
    clear(vars.normalization{:});
end

%% Comptuing the Bhattacharyya and Statistical distance
if isGetDistance == 1 && isNormalize == 1
    numFeats = size(dataCatheterNorm,2)-1;
    distBhatt = zeros(1, numFeats);
    distStat = zeros(1, numFeats);
    distType = 'mahalanobis';
    for i = 1:numFeats
        distBhatt(1, i) = GetBhattacharyyaDistance(dataCatheterNorm(:, i+1), dataTissueNorm(:, i+1)); 
        temp = pdist2(dataCatheterNorm(:, i+1), dataTissueNorm(:, i+1),...
                      distType, 'Smallest', 1); 
        distStat(1, i) = mean(temp); 
    end
    vars.distanceComputing = {'temp', 'distType', 'i', 'isGetDistances', ...
                              'numFeats', 'isGetDistance'};
    clear(vars.distanceComputing{:});
end

if exist('isGetDistance', 'var') == 1
    clear isGetDistance;
end

if exist('isNormalize', 'var') == 1
    clear isNormalize;
end
%% Get not normalized data
dataTraining = struct2cell(dataTraining.').';
% if isUseNormalization == 0
% %     dataTraining = struct2cell(dataTraining.').';
%     disp('Alexei');
% else
%     dataTraining = dataTrainingNorm;
% end
inputTrainStruct = struct('Area', dataTraining(1:end,2), ...
                          'ConvexArea', dataTraining(1:end,3), ...
                          'Perimeter', dataTraining(1:end,4), ...
                          'Eccentricity', dataTraining(1:end,5), ...
                          'Solidity', dataTraining(1:end,6), ...
                          'Extent', dataTraining(1:end,7), ...
                          'EquivDiameter', dataTraining(1:end,8), ...
                          'MaxIntensity', dataTraining(1:end,9), ...
                          'MeanIntensity', dataTraining(1:end,10), ...
                          'MinIntensity', dataTraining(1:end,11), ...
                          'Variance', dataTraining(1:end,12), ...
                          'StandardDeviation', dataTraining(1:end,13), ...
                          'Contrast', dataTraining(1:end,14), ...
                          'Correlation', dataTraining(1:end,15), ...
                          'Energy', dataTraining(1:end,16), ...
                          'Homogeneity', dataTraining(1:end,17), ...
                          'MajorAxisLength', dataTraining(1:end,18), ...
                          'MinorAxisLength', dataTraining(1:end,19));

targetTrainStruct = struct('Presence', dataTraining(1:end,1));
netTrainTargets = struct2mat(targetTrainStruct);
% targetTrainDataset = struct2dataset(targetTrainStruct); 
clear targetTrainStruct;
netTrainInputs = struct2mat(inputTrainStruct);
% inputTrainDataset = struct2dataset(inputTrainStruct); 
clear inputTrainStruct;
% clear dataTraining;
%% Check for NaNs
if CheckNaN(netTrainInputs) > 0
    msg = sprintf('There are %d NaNs in your training data!', CheckNaN(netTrainInputs));
    helpdlg(msg, 'Point Selection');
end
