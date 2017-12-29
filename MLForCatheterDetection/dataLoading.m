%% Initial state
clear all; close all; clc;
set(0, 'DefaultFigureWindowStyle', 'normal');
currentFolder = pwd;
addpath(genpath(pwd));
isGetDistance = 1;
useNormalizedVersion = 0;
isLoadSeparatedData = 0;

%% Load the data                             
dataTraining = struct('Presence', [], ...
                      'Area', [], ... 
                      'ConvexArea', [], ...
                      'Perimeter', [], ...
                      'Eccentricity', [], ...
                      'Solidity', [], ...
                      'Extent', [], ...
                      'EquivDiameter', [], ...
                      'MaxIntensity', [], ...
                      'MeanIntensity', [], ...
                      'Variance', [], ...
                      'StandardDeviation', [], ...
                      'Skewness', [], ...
                      'Kurtosis', [], ...
                      'Contrast', [], ...
                      'Correlation', [], ...
                      'Entropy', [], ...
                      'Energy', [], ...
                      'Homogeneity', [], ...
                      'MajorAxisLength', [], ... 
                      'MinorAxisLength', []);
dataTrainingNorm = dataTraining;

% Load training data

if isLoadSeparatedData == 1
    trainingFolder = strcat(currentFolder, '\Labeled data\LV Catheter 07 v.6 (20 features + separation)\Training');
else
    trainingFolder = strcat(currentFolder, '\Labeled data\LV Catheter 07 v.7 (20 features)\Training');
end
cd(trainingFolder);
s = what;
matfiles = s.mat;
for i = 1:numel(matfiles)
    tempTrain = load(char(matfiles(i)));
    tempTrainNorm = tempTrain.featuresAllNorm; 
    tempTrain = tempTrain.featuresAll;
    dataTrainingNorm = [tempTrainNorm, dataTrainingNorm]; 
    dataTraining = [tempTrain, dataTraining];
end
dataTraining(end) = [];
dataTrainingNorm(end) = [];
cd ..\..\..\
vars.loadTrainData = {'currentFolder', 'trainingFolder' , 'i', 's', ...
                      'matfiles', 'tempTrain', 'tempTrainNorm'};
clear(vars.loadTrainData{:});

%% Removing NaNs
fn = fieldnames(dataTraining);
for i = 1:numel(fn)
    dataTraining = dataTraining(~isnan([dataTraining.(fn{i})]));
    dataTrainingNorm = dataTrainingNorm(~isnan([dataTrainingNorm.(fn{i})]));
end
vars.removingNaNs = {'fn', 'i'};
clear(vars.removingNaNs{:});

%% Get normalized data
dataNorm = dataTrainingNorm; 
numFields = numel(fieldnames(dataTrainingNorm));
dataTissueNorm = zeros(1, numFields, 'double');
dataCatheterNorm = zeros(1, numFields, 'double');
dataTrainingNorm = struct2mat(dataTrainingNorm); 
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
vars.normalization = {'dataTrainingNorm', 'numFields', 'dataTemp', 'i', 'temp'};
clear(vars.normalization{:});

%% Get non-normalized data
data = dataTraining;   
numFields = numel(fieldnames(dataTraining));  
dataTissue = zeros(1, numFields, 'double');
dataCatheter = zeros(1, numFields, 'double');
dataTraining = struct2mat(dataTraining); 
for i = 1:size(dataTraining, 1)
   if dataTraining(i,1) == 1
        dataCatheter(end+1,:) = dataTraining(i,:); 
   else 
        dataTissue(end+1,:) = dataTraining(i,:);  
   end   
end
netTrainTargets = dataTraining(:,1);
netTrainInputs = dataTraining;
netTrainInputs(:,1) = [];
dataCatheter(1,:) = [];
dataTissue(1,:) = [];
vars.normalization = {'dataTraining', 'numFields', 'dataTemp', 'i', 'temp'};
clear(vars.normalization{:});

%% Comptuing Bhattacharyya and Statistical distance
if isGetDistance == 1
    
    if useNormalizedVersion == 0
        cathData = dataCatheter;
        tissueData = dataTissue;
    elseif useNormalizedVersion == 1
        cathData = dataCatheterNorm;
        tissueData = dataTissueNorm;
    end
    
    numFeats = size(cathData,2)-1;
    distBhatt = zeros(1, numFeats);
    distStat = zeros(1, numFeats);
    distType = 'euclidean';
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
if useNormalizedVersion == 0
    if CheckNaN(netTrainInputs) > 0
        msg = sprintf('There are %d NaNs in your training data!', CheckNaN(netTrainInputs));
        helpdlg(msg, 'Point Selection');
    end
elseif useNormalizedVersion == 1
    if CheckNaN(netTrainInputsNorm) > 0
        msg = sprintf('There are %d NaNs in your training data!', CheckNaN(netTrainInputs));
        helpdlg(msg, 'Point Selection');
    end
end
    
if exist('isNormalize', 'var') == 1
    clear isNormalize;
end