%   input - input data.
%   output - target data.
% x = trainInputs;
% t = trainTargets;
x = netTrainInputs';
t = netTrainTargets';
layerSize = 'small';
trainingFunction = 'BR';    % training function   
isGPU = 'no';               % train a net on GPUs
isParallel = 'no';          % train a net on a parallel pool

% Create a pool
pool = gcp('nocreate');             % train a net on a parallel pool
if strcmp(isParallel,'yes')
    parallelPool;                   % only for SVM and KNN
    disp('Parallel computatinos')
elseif strcmp(isParallel,'no')
    disp('Serial computatinos');
end

% Training Function
% For a list of all training functions type: help nntrain
% 'trainlm' is usually fastest.
% 'trainbr' takes longer but may be better for challenging problems.
% 'trainscg' uses less memory. Suitable in low memory situations.
% 'trainscg' Scaled conjugate gradient backpropagation.

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

% Create a Pattern Recognition Network
switch layerSize
    case 'small'
        hiddenLayerSize = [20, 10, 5];
    case 'mid'
        hiddenLayerSize = [40, 20, 10];
    case 'big'
        hiddenLayerSize = [180, 140, 120];
    otherwise
        disp('Choose size of network proprely')
end

net = cascadeforwardnet(hiddenLayerSize,trainFcn);

% Transfer function of the layers
numLayers = numel(hiddenLayerSize);
for i = 1:numLayers
    net.layers{i}.transferFcn = 'tansig';
end
net.layers{numLayers+1}.transferFcn = 'tansig';

% Setup Division of Data for Training, Validation, Testing
net.divideFcn = 'dividerand';  % Divide data randomly
net.divideMode = 'sample';  % Divide up every sample
net.divideParam.trainRatio = 70/100;
net.divideParam.valRatio = 15/100;
net.divideParam.testRatio = 15/100;


% Number of epochs
net.trainParam.epochs = 1000;

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
figure, plotconfusion(t,y, 'Overall')
% figure, plotroc(t,y)