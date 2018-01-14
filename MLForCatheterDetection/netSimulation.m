%% Initial state
clear all;  close all; clc;
addpath(genpath(pwd));

% Initial variables
useNormalizedData = 1;      % use normalized type of data (1) or not (0)
isLoadSeparatedData = 0;    % separated data = 1, not separated data = 0
sfsType = 'Full';           % Full, DA, SVM, KNN, FSRA, BDFS, OFS, SBFS
numFeats = 12;

% Load the data
if isLoadSeparatedData == 0
    simData = load('simData (not separated).mat');
else
    simData = load('simData (separated).mat');
end

if useNormalizedData == 1
    inputs = simData.netTrainInputsNorm';
    targets = simData.netTrainTargetsNorm';
elseif useNormalizedData == 0
    inputs = simData.netTrainInputs';
    targets = simData.netTrainTargets';
end

% Choose SFS model
switch sfsType
    case 'Full'
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
[simInputs, numFeatures] = GetDataUsingModel(inputs', isLoadSeparatedData, modelRange);
simInputs = simInputs';
simTargets = targets';
%%
% Load net object HERE
%%
% Network application on the independent data 
[numRows, numCols] = size(simInputs);
simOutputs = zeros(numCols,1);
tic;
for i = 1:numCols
    simOutputs(i,:) = round(net(simInputs(:,i)));
end
toc;

[~, confMatrix] = confusion(simTargets', simOutputs');
confMatrix = confMatrix';

numRightTissueCases = confMatrix(1,1); 
numWrongTissueCases = confMatrix(2,1);
numRightCathCases = confMatrix(2,2); 
numWrongCathCases = confMatrix(1,2);
numTissueCases = numRightTissueCases + numWrongTissueCases;
numCathCases = numRightCathCases + numWrongCathCases;

percentTissueAccuracy = numRightTissueCases / numTissueCases;
percentTissueErrors = 1 - percentTissueAccuracy;
percentCathAccuracy = numRightCathCases / numCathCases;
percentCathErrors = 1 - percentCathAccuracy;

fprintf('\n*************************************************************************')
fprintf('\n \t \t Catheter classification rate: %.3f (%.d right cases out of %.d) \n', ...
         percentCathAccuracy, numRightCathCases, numCathCases)
fprintf('\n \t \t Catheter misclassification rate: %.3f (%.d wrong cases out of %.d) \n', ...
        percentCathErrors, numWrongCathCases, numCathCases)
fprintf('\n \t \t Tissue classification rate: %.3f (%.d right cases out of %.d) \n', ...
         percentTissueAccuracy, numRightTissueCases, numTissueCases)
fprintf('\n \t \t Tissue misclassification rate: %.3f (%.d wrong cases out of %.d) \n', ...
        percentTissueErrors, numWrongTissueCases, numTissueCases)   
fprintf('************************************************************************* \n')


