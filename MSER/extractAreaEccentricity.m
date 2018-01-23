%% Initial state of the programm
clear all; close all; clc;
set(0, 'DefaultFigureWindowStyle', 'normal');
currentFolder = pwd;
addpath(genpath(pwd));

%% Global variables
isVisual = 0;
isFill = 0;
nTimeframe = 9; %9
sliceRange = 108:110; % 35:55
scrSz = get(0, 'Screensize');
featuresAll = struct([]);

%% Reading the data
tic;
filename = 'LV Catheter 07.nrrd'; % delete after getting  features
[X, meta] = nrrdread(filename);
Y = double(X);
sz = sscanf(meta.sizes, '%d');
nDims = sscanf(meta.dimension, '%d');
I = squeeze(X(:,:,:,nTimeframe));
% implay(I);
toc;

%% Binarization
for nSlice = sliceRange % comment if you need to check 1 slice
    img = I(:,:,nSlice);
    % level = threshTool(img)/255;
    [level,EM] = graythresh(img);
    BW = imbinarize(img, level);
    
    %% Filling holes
    if isFill == 1
        BWfill = imfill(BW, 'holes');
    elseif isFill == 0
        BWfill = BW;
    end
    CCfill = bwconncomp(BWfill);
    Lfill = labelmatrix(CCfill);
    numFill = CCfill.NumObjects;
        
    if isVisual == 1
        colorLabelIni = label2rgb(Lfill, 'parula', 'k', 'shuffle');
        str1 = sprintf('Region labeling');
        str2 = sprintf('Objects found: %d', numFill);
        imshow(colorLabelIni, 'InitialMagnification', 'fit'); addTitle({str1; str2});
        set(gcf, 'Position', scrSz, 'Color', 'w');
    end
    vars.fillingHoles = {'str1', 'str2'};
    clear(vars.fillingHoles{:});
    
    featuresFill = regionprops(Lfill, {'Area', 'Eccentricity', 'Extrema'});

        if isVisual == 0
        hFig = figure;
%         imshow(BWfill, 'InitialMagnification', 'fit');
        imshowpair(BWfill, img, 'montage');
        str1 = sprintf("%d slice", nSlice);
        str2 = 'Enumeration of regions';
        addTitle({str1, str2});
        hold on
        for count = 1:numFill
            posX = featuresFill(count).Extrema(1,1) - 4;
            posY = featuresFill(count).Extrema(1,2) - 7;
            str = sprintf("%d", count);
            text(posX, posY, char(str), 'FontSize', 18, 'FontName', 'Times New Roman', 'Color', 'g');
        end
%         set(gcf, 'Position', scrSz/2, 'Color', 'w', 'name', str1, 'numbertitle', 'off');
        set(gcf, 'Position', [scrSz(3)/2, scrSz(2), scrSz(3)/2, scrSz(4)],...
            'Color', 'w', 'name', str1, 'numbertitle', 'off');
        hold off
    end
    close(hFig);
    disp(str1);
   
end % comment if you need to check 1 slice
