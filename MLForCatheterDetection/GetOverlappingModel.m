%% Initial state
clear all;  close all; clc;
addpath(genpath(pwd));
epsilon = 0.05;

% Initial variables
isLoadSeparatedData = 0;    % separated data = 1, not separated data = 0

% Load the data
if isLoadSeparatedData == 0
    data = load('data (not separ.).mat');
else
    data = load('data (sep.).mat');
end

dataCatheter = data.dataCatheter;
dataTissue = data.dataTissue;
dataCatheter(:,1) = [];
dataTissue(:,1) = [];

numFeats = size(dataCatheter,2);
overlapping = zeros(1, numFeats);

for i = 1:numFeats
    overlapping(1,i) = GetOverlapRate(dataCatheter(:,i), dataTissue(:,i), epsilon)/100;
end
