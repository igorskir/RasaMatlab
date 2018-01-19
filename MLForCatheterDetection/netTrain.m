%% DO NOT NORMOLIZE TARGRETS
%% Initial state
clear all;  close all; clc;
addpath(genpath(pwd));

% Initial variables
isVisual = 0;
useNormalizedData = 1;      % use normalized type of data (1) or not (0)
sfsType = 'Full';           % Full, DA, SVM, KNN, FSRA, BDFS, OFS, SBFS
numFeats = 6;               % 20 (Full), 12 and 6    
netType = 'feed-forward';   % 'feed-forward', 'cascade', 'recurrent'
netSize = 'mid';            % small, mid, big
trainingFunction = 'BR';    % training function  
isGPU = 'no';               % train a net on GPUs
isParallel = 'no';          % train a net on a parallel pool
isLoadSeparatedData = 0;    % separated data = 1, not separated data = 0

% Load the data
if isLoadSeparatedData == 0
    x = load('inputs (not separated).mat');
    t = load('targets (not separated).mat');
else
    x = load('inputs (separated).mat');
    t = load('targets (separated).mat');
end

if useNormalizedData == 1
    x = x.netTrainInputsNorm';
    t = t.netTrainTargetsNorm';
elseif useNormalizedData == 0
    x = x.netTrainInputs';
    t = t.netTrainTargets';
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
[x, numFeatures] = GetDataUsingModel(x', isLoadSeparatedData, modelRange);
x = x';
% Use of the particular model
% model = zeros(1,20);
% model(1,1) = 1;
% x = GetDataUsingModel(x', model)';

% Create a pool
pool = gcp('nocreate');             
if strcmp(isParallel,'yes')
    parallelPool;                   % only for SVM and KNN
    disp('Parallel computations')
elseif strcmp(isParallel,'no')
    disp('Serial computations');
end

% Choose a training fucntion
switch trainingFunction
    case 'LM'
        trainFcn = 'trainlm';
    case 'BR'
        trainFcn = 'trainbr';
    case 'BFG'
        trainFcn = 'trainbfg';
    case 'RP'
        trainFcn = 'trainrp';
    case 'SCG'
        trainFcn = 'trainscg';
    case 'CGB'
        trainFcn = 'traincgb';
    case 'CGF'
        trainFcn = 'traincgf';
    case 'CGP'
        trainFcn = 'traincgp';
    case 'OSS'
        trainFcn = 'trainoss';
    case 'GDX'
        trainFcn = 'traingdx';
    otherwise
        disp('Choose the training methods proprely')
end

% Choose a layer size
switch netSize
    case 'small'
        hiddenLayerSize = [20, 10, 5];
    case 'mid'
        hiddenLayerSize = [40, 20, 10];
    case 'big'
        hiddenLayerSize = [60, 30, 15];
    otherwise
        disp('Choose the size of network proprely')
end

% Performance function
if strcmp (trainFcn, 'trainscg')
   performFcn = 'crossentropy';
else
   performFcn = 'mse'; 
end

% Create a network object
numLayers = numel(hiddenLayerSize);
switch netType
    case 'feed-forward'
        net = patternnet(hiddenLayerSize, trainFcn, performFcn);
%         net.layers{numLayers+1}.transferFcn = 'logsig';   
        net.layers{numLayers+1}.transferFcn = 'tansig';      
    case 'cascade'
        net = cascadeforwardnet(hiddenLayerSize,trainFcn);
        net.layers{numLayers+1}.transferFcn = 'tansig';        
    case 'recurrent'
        layerDelays = 1:2;
        net = layrecnet(layerDelays,hiddenLayerSize,trainFcn);
        net.layers{numLayers+1}.transferFcn = 'tansig';
    otherwise
        disp('Choose the network type proprely');
end

% Transfer function of hidden layers
for i = 1:numLayers
    net.layers{i}.transferFcn = 'tansig';
end

% Setup Division of Data for Training, Validation, Testing
net.divideFcn = 'dividerand';  % Divide data randomly
net.divideMode = 'sample';  % Divide up every sample
net.divideParam.trainRatio = 70/100;
net.divideParam.valRatio = 15/100;
net.divideParam.testRatio = 15/100;
net.trainParam.max_fail = 10;

% Number of epochs
net.trainParam.epochs = 1000;

% Train a network
if (strcmp(isParallel,'yes') && strcmp(isGPU, 'yes'))
    error('Choose only one mode of training')
end

if strcmp(isParallel, 'yes')
    [net,tr] = train(net,x,t, 'useParallel', 'yes', 'showResources', 'yes');
    delete(gcp('nocreate'))
elseif strcmp(isGPU,'yes')
    [net,tr] = train(net,x,t, 'useGPU', 'yes');
else
    [net,tr] = train(net,x,t);
end

% Test the Network
y = net(x);
e = gsubtract(t,y);
performance = perform(net,t,y)
tind = vec2ind(t);
yind = vec2ind(y);
percentErrors = sum(tind ~= yind)/numel(tind);

% View the Network
if isVisual == 1
    view(net)
end

% Plots

% Training Confusion Plot Variables
yTrn = net(x(:,tr.trainInd));
tTrn = t(:,tr.trainInd);

% Validation Confusion Plot Variables
yVal = net(x(:,tr.valInd));
tVal = t(:,tr.valInd);

% Test Confusion Plot Variables
yTst = net(x(:,tr.testInd));
tTst = t(:,tr.testInd);

% Overall Confusion Plot Variables
yAll = net(x);
tAll = t;

% Plot Confusion
if isVisual == 1
    figure, plotconfusion(tTrn, yTrn, 'Training', ...
                          tVal, yVal, 'Validation', ...
                          tTst, yTst, 'Test', ...
                          tAll, yAll, 'Overall')
    figure, plotconfusion(t,y)
    % Uncomment these lines to enable various plots.
    % figure, plotperform(tr)
    % figure, plottrainstate(tr)
    % figure, ploterrhist(e)
    % figure, plotroc(t,y)
end

% Get the results
[~, confMatrix] = confusion(tAll,yAll);
confMatrix = confMatrix';

accuracyVals = zeros(1,10);
accuracyVals(1,1) = confMatrix(1,1)/sum(confMatrix(:,1));
accuracyVals(1,2) = confMatrix(2,2)/sum(confMatrix(:,2));
accuracyVals(1,3) = 1 - 1/((confMatrix(1,1) + confMatrix(2,2))/(confMatrix(2,1) + confMatrix(1,2)));
accuracyVals(1,4) = 0;
accuracyVals(1,5) = 1 - accuracyVals(1,1);
accuracyVals(1,6) = 1 - accuracyVals(1,2);
accuracyVals(1,7) = 1 - accuracyVals(1,3);
accuracyVals(1,8) = 0;
accuracyVals(1,9) = performance;
accuracyVals(1,10) = tr.time(end);
stopVal = tr.stop;
netResult = {accuracyVals stopVal};

fprintf('Catheter Classification Rate: %.2f%%\n', 100*accuracyVals(1,2));
fprintf('Catheter Misclassification Rate: %.2f%%\n', 100*accuracyVals(1,6));

% Save net-file
layersName = [num2str(hiddenLayerSize(1,1)), ' ', num2str(hiddenLayerSize(1,2)), ' ', num2str(hiddenLayerSize(1,3))];
if strcmp(sfsType, 'Full')
    netFilename = [netType, ' ', num2str(numFeatures), '_features [', layersName, ']'];
else
    netFilename = [netType, ' ', num2str(numFeatures), '_features [', layersName, '] ', sfsType];
end
cd('Net tests')
save(netFilename, 'net', 'tr', 'x', 'y', 't');
cd ..\