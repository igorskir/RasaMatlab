%% DO NOT NORMOLIZE TARGRETS
%% Initial state
clear all; close all; clc;

% Initial variables
isNormalize = 1;
numFeats = 18; % 5, 9, 11, 18

% Load the data
if isNormalize == 1 
    load(strcat(num2str(numFeats), '_features_norm.mat'));
    load('targets');
elseif isNormalize == 0
    load(strcat(num2str(numFeats), '_features.mat'));
    load('targets');
end

x = netSelectInputs';
%x = netTrainInputs';
t = netTrainTargets';

netType = 'feed-forward'; % 'feed-forward', 'cascade', 'recurrent'
netSize = 'small'; % small, mid , big
trainingFunction = 'BR';    % training function  
isGPU = 'no';               % train a net on GPUs
isParallel = 'no';          % train a net on a parallel pool

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
view(net)

% Plots
% Uncomment these lines to enable various plots.
% figure, plotperform(tr)
% figure, plottrainstate(tr)
% figure, ploterrhist(e)

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
figure, plotconfusion(tTrn, yTrn, 'Training', ...
                      tVal, yVal, 'Validation', ...
                      tTst, yTst, 'Test', ...
                      tAll, yAll, 'Overall')
figure, plotconfusion(t,y)
% figure, plotroc(t,y)

