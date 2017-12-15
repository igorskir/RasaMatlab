%% Initial state of the programm
clear all; close all; clc;
set(0, 'DefaultFigureWindowStyle', 'normal');
currentFolder = pwd;
addpath(genpath(pwd));

%% Global variables
isVisual = 1;
isSeparate = 1;
isTestThreshLevel = 1;
openArea = 15;
method = 'euclidean';
nTimeframe = 13; %9
k = 1.1;
ax = 'short'; % 'long1', 'long2'
   
%XLS reading
cathDataFile = 'LV Catheter 07.xlsx';
sheet = 'Sheet1';
switch ax
    case 'short'
        xlRange = 'C5:S6';
    case 'long1'
        xlRange = 'C15:S16';
    case 'long2'
        xlRange = 'C25:S26';
    otherwise
        disp('Wrong value')
end

cathData = xlsread(cathDataFile,sheet,xlRange);
minSlice = cathData(1,nTimeframe);
maxSlice = cathData(2,nTimeframe);
sliceRange = minSlice:maxSlice;
scrSz = get(0, 'Screensize');

%% Reading the data
filename = 'LV Catheter 07.nrrd';
[X, meta] = nrrdread(filename);
sz = sscanf(meta.sizes, '%d');
nDims = sscanf(meta.dimension, '%d');
I = squeeze(X(:,:,:,nTimeframe));
%% Estimate threshold levels 
if isTestThreshLevel == 1
testThreshLevels = zeros(numel(sliceRange),2);
    for nSlice = sliceRange
%         img = GetImage(I, nSlice, ax);
        img = procImage;
        levelAuto = threshTool(img);
        levelMan = graythresh(img);
        testThreshLevels(nSlice, 1) = levelAuto;
        testThreshLevels(nSlice, 2) = levelMan; 
    end
end
%% Compare BWs with 2 different levels
for nSlice = sliceRange
%     img = GetImage(I, nSlice, ax);
    img = procImage;
    [level,EM] = graythresh(img);
    BW1 = imbinarize(img, level);
    BW2 = imbinarize(img, k*level);
    BW1 = bwareaopen(BW1, openArea);
    BW2 = bwareaopen(BW2, openArea);
    if isSeparate == 1
        BW1 = GetSeparatedRegions(BW1, method, ax);
        BW2 = GetSeparatedRegions(BW2, method, ax);
    end
        if isVisual == 1
            str1 = sprintf('Binarized image: %d', nSlice);
            str2 = sprintf('Thresholding levels: %.3f and %.3f (%d and %d out of 255)', ...
                            level, k*level, uint8(level*255), uint8(k*level*255));
%             imshowpair(BW1, BW2, 'montage');
            imshowpair(BW2, img, 'montage');
            AddTitle({str1; str2});
            set(gcf, 'Position', [scrSz(3), 0, scrSz(3), scrSz(4)],...
            'Color', 'w', 'name', str1, 'numbertitle', 'off');
%             set(gcf, 'Position', scrSz, 'Color', 'w');
        end
    
end