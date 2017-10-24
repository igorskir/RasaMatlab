%% Initial state of the programm
clear all; close all; clc;
set(0, 'DefaultFigureWindowStyle', 'normal');
currentFolder = pwd;
addpath(genpath(pwd));

%% Keep only ones in "Presence" column
featuresOnes = struct('Presence', [], 'Area', [], ... 
                      'Eccentricity', [], 'EquivDiameter', [], ...
                      'MaxIntensity', [], 'MeanIntensity', [], ...
                      'MinIntensity', [], 'Variance', [], ...
                      'StandardDeviation' , [], 'Contrast', [], ...
                      'Correlation', [], 'Energy', [], ...
                      'Homogeneity', [], 'MajorAxisLength', [], ... 
                      'MinorAxisLength', [], 'Orientation', []);

temp = featuresAll;
numStrings = numel(temp);
for i = 1:numStrings %numStrings
    if temp(i).Presence == 1
       featuresOnes = [temp(i), featuresOnes]; 
    end
end
featuresOnes(end) = [];