%% DO NOT NORMOLIZE TARGRETS
%% Initial state
clear all; close all; %clc;
set(0, 'DefaultFigureWindowStyle', 'normal');
currentFolder = pwd;
addpath(genpath(pwd));
isNormalize = 1;
netType = 'feed-forward'; % 'feed-forward', 'cascade', 'recurrent'
netSize = 'big'; %'big', 'mid', 'small'
%% Load simulation data
dataSimulation = struct('Presence', [], ...
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

simulationFolder = strcat(currentFolder, '\Labeled data\LV Catheter 07 v.3 (18 features)\Simulation');
cd(simulationFolder);
s = what;
matfiles = s.mat;
for i = 1:numel(matfiles)
    tempSim = load(char(matfiles(i)));
    tempSim = tempSim.featuresAll;
    dataSimulation = [tempSim, dataSimulation];
end
dataSimulation(end) = [];
cd ..\..\..\

vars.loadSimData = {'currentFolder','simulationFolder', 'i', 's', ...
                    'matfiles', 'tempSim'};
clear(vars.loadSimData{:});

%% Removing NaNs
fn = fieldnames(dataSimulation);
for i = 1:numel(fn) 
    dataSimulation = dataSimulation(~isnan([dataSimulation.(fn{i})]));
end
vars.removingNaNs = {'fn', 'i'};
clear(vars.removingNaNs{:});

%%  Prepare simulation data
dataSimulation = struct2cell(dataSimulation.').';
inputSimStruct = struct('Area', dataSimulation(1:end,2), ...
                        'ConvexArea', dataSimulation(1:end,3), ...
                        'Perimeter', dataSimulation(1:end,4), ...
                        'Eccentricity', dataSimulation(1:end,5), ...
                        'Solidity', dataSimulation(1:end,6), ...
                        'Extent', dataSimulation(1:end,7), ...
                        'EquivDiameter', dataSimulation(1:end,8), ...
                        'MaxIntensity', dataSimulation(1:end,9), ...
                        'MeanIntensity', dataSimulation(1:end,10), ...
                        'MinIntensity', dataSimulation(1:end,11), ...
                        'Variance', dataSimulation(1:end,12), ...
                        'StandardDeviation', dataSimulation(1:end,13), ...
                        'Contrast', dataSimulation(1:end,14), ...
                        'Correlation', dataSimulation(1:end,15), ...
                        'Energy', dataSimulation(1:end,16), ...
                        'Homogeneity', dataSimulation(1:end,17), ...
                        'MajorAxisLength', dataSimulation(1:end,18), ...
                        'MinorAxisLength', dataSimulation(1:end,19));

targetSimStruct = struct('Presence', dataSimulation(1:end,1));
netSimInputs = struct2mat(inputSimStruct); clear inputSimStruct;
netSimTargets = struct2mat(targetSimStruct); clear targetSimStruct;

if isNormalize == 1
    netSimInputs = robustNormalization(netSimInputs, 'quantile', 1);
%     netSimTargets = robustNormalization(netSimTargets, 'quantile');
end

vars.loadSimData = {'inputSimStruct;', 'targetSimStruct', 'isNormalize'};
clear(vars.loadSimData{:});

%% Check for NaNs
if checkNaN(netSimInputs) > 0
    msg = sprintf('There are %d NaNs in your training data!', checkNaN(netSimInputs));
    helpdlg(msg, 'Point Selection');
end

%% Load a net object
netFolder = 'D:\Clouds\Google Drive\RASA\Matlab\12. Catheter ML\Net tests\Normolized data\';
switch netType  
    case 'feed-forward'
        netTypeFolder = '1. Feed-forward backprop';       
    case 'cascade'
        netTypeFolder = '2. Cascade-forward';        
    case 'recurrent'
        netTypeFolder = '3. Reccurent'; 
    otherwise
        disp('Choose the network type proprely');
end

switch netSize
    case 'small'
        netSizeFolder = '[20 10 5]';
    case 'mid'
        netSizeFolder = '[40 20 10]';
    case 'big'
        netSizeFolder = '[80 40 20]';
    otherwise
        disp('Choose the layer layout properly')
end
load(fullfile(netFolder, netTypeFolder,netSizeFolder,'net.mat'));
netTypeFolder(1:3) = [];
%% Results of net application
[numRows, numCols] = size(netSimInputs);
netSimInputs = netSimInputs';
netSimResults = zeros(numRows,1);
tic;
for i = 1:numRows
    netSimResults(i,:) = round(net(netSimInputs(:,i)));
end
toc;

%% Perfomance
numClasCases = numRows - sum(netSimResults ~= netSimTargets);
netSimAccuracy = 100*sum(netSimResults == netSimTargets)/numel(netSimTargets);
numMisclasCases  = numRows - sum(netSimResults == netSimTargets);
netSimError  = 100*sum(netSimResults ~= netSimTargets)/numel(netSimTargets);

fprintf('****************************************************************')
fprintf('\n \n \t \t %s %s network \n', netTypeFolder, netSizeFolder)
fprintf('\n \t \t Classification rate: %.2f (%.d right cases) \n', ...
         netSimAccuracy, numClasCases)
fprintf('\n \t \t Misclassification rate: %.2f (%.d wrong cases) \n \n', ...
        netSimError, numMisclasCases)
fprintf('**************************************************************** \n')

