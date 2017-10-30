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
trainingFolder = strcat(currentFolder, '\Labeled data\LV Catheter 07 v.3 (18 features)\Training');
% trainingFolder = strcat(currentFolder, '\Labeled data\LV Catheter 07 v.4 (18 features)\Training');
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
    temp = struct2mat(dataTraining);
    dataTrainingNorm = RobustNormalization(temp, 'quantile', 0); % quartile, std, mean, linear
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
    vars.normalization = {'numFields', 'dataTemp', 'i', 'temp'};
    clear(vars.normalization{:});
end

%% Get not normalized data
if isNormalize == 0
    temp = struct2mat(dataTraining);
    numFields = numel(fieldnames(dataTraining));
    dataTissue = zeros(1, numFields, 'double');
    dataCatheter = zeros(1, numFields, 'double');
    for i = 1:size(temp, 1)
       if temp(i,1) == 1
            dataCatheter(end+1,:) = temp(i,:); 
       else 
            dataTissue(end+1,:) = temp(i,:);  
       end   
    end
    netTrainTargets = temp(:,1);
    netTrainInputs = temp;
    netTrainInputs(:,1) = [];
    dataCatheter(1,:) = [];
    dataTissue(1,:) = [];
    vars.normalization = {'numFields', 'dataTemp', 'i', 'temp'};
    clear(vars.normalization{:});
end

%% Comptuing Bhattacharyya and Statistical distance
if isGetDistance == 1
    
    if isNormalize == 0
        cathData = dataCatheter;
        tissueData = dataTissue;
    elseif isNormalize == 1
        cathData = dataCatheterNorm;
        tissueData = dataTissueNorm;
    end
    
    numFeats = size(cathData,2)-1;
    distBhatt = zeros(1, numFeats);
    distStat = zeros(1, numFeats);
    distType = 'mahalanobis';
    for i = 1:numFeats
        distBhatt(1, i) = GetBhattacharyyaDistance(cathData(:, i+1), tissueData(:, i+1)); 
        temp = pdist2(cathData(:, i+1), tissueData(:, i+1),...
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

%% Check for NaNs

if isNormalize == 0
    if CheckNaN(netTrainInputs) > 0
        msg = sprintf('There are %d NaNs in your training data!', CheckNaN(netTrainInputs));
        helpdlg(msg, 'Point Selection');
    end
elseif isNormalize == 1
    if CheckNaN(netTrainInputsNorm) > 0
        msg = sprintf('There are %d NaNs in your training data!', CheckNaN(netTrainInputs));
        helpdlg(msg, 'Point Selection');
    end
end
    
if exist('isNormalize', 'var') == 1
    clear isNormalize;
end