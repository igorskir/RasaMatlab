%% Initial state
clear all;  close all; clc;
addpath(genpath(pwd));

% Initial variables
isVisual = 0;
useNormalizedData = 1;      % use normalized type of data (1) or not (0)
isLoadSeparatedData = 0;    % separated data = 1, not separated data = 0
sfsType = 'SBFS';
numFeats = 12;

% Load the data
if isLoadSeparatedData == 0
    simData = load('simData (not separated).mat');
else
    simData = load('simData (separated).mat');
end

if useNormalizedData == 1
    simInputs = simData.netTrainInputsNorm';
    simTargets = simData.netTrainTargetsNorm';
elseif useNormalizedData == 0
    simInputs = simData.netTrainInputs';
    simTargets = simData.netTrainTargets';
end

% Choose SFS model
switch sfsType
    case 'FULL'
         modelRange = 'D5:W5';
    case 'DA'
        if numFeats == 12
            modelRange = 'D6:W6';
        elseif numFeats == 6
            modelRange = 'D17:W17';
        end
    case 'SVM'
        if numFeats == 12
            modelRange = 'D7:W7';
        elseif numFeats == 6
            modelRange = 'D18:W18';
        end
    case 'KNN'
        if numFeats == 12
            modelRange = 'D8:W8';
        elseif numFeats == 6
            modelRange = 'D19:W19'; 
        end
    case 'FSRA'
        if numFeats == 12
            modelRange = 'D9:W9';
        elseif numFeats == 6
            modelRange = 'D20:W20';  
        end
    case 'BDFS'
        if numFeats == 12
            modelRange = 'D10:W10';
        elseif numFeats == 6
            modelRange = 'D21:W21'; 
        end
    case 'OFS'
        if numFeats == 12
            modelRange = 'D11:W11';
        elseif numFeats == 6
            modelRange = 'D22:W22';  
        end
    case 'SBFS'
        if numFeats == 12
            modelRange = 'D12:W12';
        elseif numFeats == 6
            modelRange = 'D23:W23';  
        end
end
        
% Get the data based on the used model
[modelInputs, numFeatures] = GetDataUsingModel(simInputs', isLoadSeparatedData, modelRange);
modelInputs = modelInputs';
modelTargets = simTargets;

