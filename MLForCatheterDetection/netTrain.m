%% DO NOT NORMOLIZE TARGRETS
%% Initial state
% clear all;  close all; clc;
addpath(genpath(pwd));

% Initial variables
isVisual = 0;
useNormalizedData = 1;   % use normalized type of data (1) or not (0)
netType = 'feed-forward';   % 'feed-forward', 'cascade', 'recurrent'
netSize = 'small';          % small, mid , big
trainingFunction = 'BR';    % training function  
isGPU = 'no';               % train a net on GPUs
isParallel = 'no';          % train a net on a parallel pool

% Load the data
x = load('inputs.mat');
t = load('targets.mat');
if useNormalizedData == 1
    x = x.netTrainInputsNorm';
    t = t.netTrainTargetsNorm';
elseif useNormalizedData == 0
    x = x.netTrainInputs';
    t = t.netTrainTargets';
end

% Get the data based on the used model
% x = GetDataUsingModel(x', 'D115:W115')';
% Use of the particular model
% model = zeros(1,20);
% model(1,1) = 1;
x = GetDataUsingModel(x', model)';

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
        hiddenLayerSize = [80, 40, 20];
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
net.trainParam.max_fail = 1000;

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

% Get classification rate
[~, confMatrix] = confusion(tAll,yAll);
confMatrix = confMatrix';
cathClassRate = confMatrix(2,2)/sum(confMatrix(:,2));
fprintf('Catheter Classification Rate: %.2f%%\n', 100*cathClassRate);
fprintf('Percentage Incorrect Classification : %.2f%%\n', 100*(1 - cathClassRate));


% EXTRA
% 
% rngBin12 = 'D111:W111'; % 62.3%
% rngDA12 = 'D44:W44'; % 73.4%
% rngSVM12 = 'D45:W45'; % 76.9%
% rngKNN12 = 'D46:W46'; % 74.6%
% rngFSRA12 = 'D47:W47'; % 75.5%
% 
% rngBin6 = 'D111:W111'; % 64.4%
% rngDA6 = 'D51:W51'; % 59.4%
% rngSVM6 = 'D52:W52'; % 60.4%
% rngKNN6 = 'D53:W53'; % 68.6%
% rngFSRA6 = 'D54:W54'; % 45.5%
