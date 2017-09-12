% Read PET data from dicom directory
fileFolder = fullfile(pwd, 'Gorohov');
functionPath = fullfile(pwd, 'Medical Image Reader and Viewer');
addpath(functionPath);
CT = readImages(fileFolder);
% PET volume is already aligned to scanner coordinate axes

%% View
tic; VolumeViewer3D(CT); toc