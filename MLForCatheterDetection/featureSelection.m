%% Feature selection
%% Initial state
clear all; % close all; clc;
addpath(genpath(pwd));
clear currentFolder;

useNormalizedData = 0; % 1 - yes, 0 - no (as IP said default value is euqal to 0)
funType = 1; % 1 - dicriminant, 2 - svm, 3 - knn, 4 - fsra
cvType = 1; % 1 - k-fold, 2 - Holdout
featType = 1; % 1 - auto, 2 - 6 featues, 3 - 12 features
options = statset('display', 'iter', 'MaxIter', 1000);
direction = 'forward'; %backward of forward
isLoadSeparatedData = 0;    % separated data = 1, not separated data = 0

% Load the data
if isLoadSeparatedData == 0
    x = load('inputs (not separated).mat');
    y = load('targets (not separated).mat');
else
    x = load('inputs (separated).mat');
    y = load('targets (separated).mat');
end

if useNormalizedData == 1
    x = x.netTrainInputsNorm;
    y = y.netTrainTargetsNorm;
elseif useNormalizedData == 0
    x = x.netTrainInputs;
    y = y.netTrainTargets;
end

% Choose number of features
switch featType
    case 1
        featTypeStr = 'Auto';
        numFeatures = [];
    case 2
        featTypeStr = '6';
        numFeatures = 6;
    case 3
        featTypeStr = '12';
        numFeatures = 12;        
end

tic;
switch cvType
    case 1
        crossValidation = cvpartition(y, 'k', 10);
        cvTypeStr = '10-fold CV';
    case 2
        crossValidation = cvpartition(y, 'Holdout', 0.2);
        cvTypeStr = '20%-holdout CV';
%         xTrain = x(crossValidation.training,:);
%         yTrain = y(crossValidation.training);
%         xTest = x(crossValidation.test,:);
%         yTest = y(crossValidation.test);
end

switch funType
    case 1
        fun = @(xTrain, yTrain, xTest, yTest) sum(yTest ~= classify(xTest, ... 
                xTrain, yTrain, 'mahalanobis'));
        funTypeStr = 'Discriminant analysis';
    case 2 
        fun = @(xTrain, yTrain, xTest, yTest) sum(yTest ~= predict(fitcsvm(xTrain, ... 
                yTrain, 'Standardize', true, 'KernelFunction', 'RBF', ...
                'KernelScale', 'auto'), xTest));
        funTypeStr = 'SVM';
    case 3
        fun = @(xTrain, yTrain, xTest, yTest) sum(yTest ~= predict(fitcknn(xTrain, ... 
                yTrain, 'Distance', 'seuclidean'), xTest));
        funTypeStr = 'KNN';
    case 4
        [~,~,~, fsraModel,~,~,fsraHistory] = stepwisefit(x, y, ...
                                                         'display', 'on');
        if featType ~= 1
            if numFeatures <= size(fsraHistory.in, 1)
                fsraModel = fsraHistory.in(numFeatures, :);
            else
                fsraModel = fsraHistory.in(end, :);
            end
        end
        funTypeStr = 'FSRA';
        
% %% Forward Stepwise Regression Algorithm
% % Forward stepwise model selection algorithm: 
% % Variables are sequentially added to the active set of variables. 
% % The procedure does not involve any tests of statistical significance 
% % of the potential covariates; instead, it produces a ranking that 
% % corresponds to the order of variables as they enter the active set.
% % The function can be provided with a dataset of size (K+1)xN 
% % (N observations, K predictors, one explained variable).
% 
% FSRARanking = fsra(netTrainTargets,netTrainInputs);

    otherwise
        msg = sprintf('Select the function type properly!');
        helpdlg(msg, 'Point Selection');
end

% Procedures for DA, SVM and KNN
if funType ~= 4
    [SFSModel, SFShistory] = sequentialfs(fun, x, y, 'cv', crossValidation, ...
                                  'options', options, 'direction', direction, ...
                                  'nfeatures', numFeatures);
    processingTime = toc;
    accSFS = 1 - SFShistory.Crit(1:end);
    numFeats = numel(accSFS);
    for i = 1:numFeats
        fprintf('Step %d, Accuracy = %.2f \n', i, accSFS(i)*100);
        pause(0.2);
    end
    fprintf('%s with %s (%s features) \n', ...
        funTypeStr, cvTypeStr, featTypeStr);
    fprintf('Processing time: %f \n', processingTime);
    fprintf('Final number of features is %d \n', numFeats);
    SFSModel = uint8(SFSModel);
end


% Procedures for FSRA
if funType == 4
    processingTime = toc;
    numFeats = sum(fsraModel);
    fprintf('%s with default cross-validation (%s features) \n', ...
            funTypeStr, featTypeStr);
    fprintf('Processing time: %f \n', processingTime);
    fprintf('Final number of features is %d \n', numFeats);
    fsraModel = uint8(fsraModel);
end

%% Extract relevant input feature array (automatic SFS)
[numRows, numCols] = size(x);
netSelectInputs = zeros(numRows, 1);
for i = 1:numCols
    if funType ~= 4
        if SFSModel(1, i) == 1
            netSelectInputs = [netSelectInputs, x(1:end, i)];
        end
    elseif funType == 4
        if fsraModel(1, i) == 1
            netSelectInputs = [netSelectInputs, x(1:end, i)];
        end 
    end
end
netSelectInputs(:,1) = [];

%% Extract relevant input feature array (manual SFS)
filename = 'D:\RASA Lab\MLForCatheterDetection\!Calculations\Feature engineering.xlsm';
sheet = 6;
useNormalizedData = 1;
rng = 'all'; % all - 18, big - 11, mid - 9, small - 5
switch rng
    case 'all'
        xlRange = 'D111:W111';
    case 'big'
        xlRange = 'D111:U111';  
    case 'mid'
        xlRange = 'D112:U112';
    case 'small'
        xlRange = 'D113:U113';
end

% Load the data
x = load('inputs.mat');
y = load('targets.mat');
if useNormalizedData == 1
    x = x.netTrainInputsNorm;
    y = y.netTrainTargetsNorm;
elseif useNormalizedData == 0
    x = x.netTrainInputs;
    y = y.netTrainTargets;
end

model = xlsread(filename, sheet, xlRange);
[numRows, numCols] = size(x);
netSelectInputs = zeros(numRows, 1);
for i = 1:numCols
        if model(1, i) == 1
            netSelectInputs = [netSelectInputs, x(1:end, i)];
        end
end
netSelectInputs(:,1) = [];
vars.featureSelection = {'filename', 'sheet', 'rng', 'model', 'numRows',...
                         'numCols','i'};
clear(vars.featureSelection{:});
