%% Initial state of the programm
clear all; close all; clc;
set(0, 'DefaultFigureWindowStyle', 'normal');
currentFolder = pwd;
addpath(genpath(pwd));

%% Global variables
isSeparate = 1; % perform separation of catheter and pigtail (1) or not
isLabel = 1; % labeling (1) or without labeling (0)
isNormalize = 1; % normalization (1) or no normaliztion (0)
isVisual = 1; % visualize (1) or not(0)
isFill = 0; % filling the holes (1) or not (0)
openArea = 15;
method = 'euclidean';

nTimeframe = 9; % 9 tf and 49/92 slice
ax = 'short'; % 'long1', 'long2'
    
%XLS reading
cathDataFile = 'LV Catheter 07.xlsx';
sheet = 'Sheet1';
switch ax
    case 'short'
        xlRange = 'C7:S8';
    case 'long1'
        xlRange = 'C17:S18';
    case 'long2'
        xlRange = 'C27:S28';
    otherwise
        disp('Wrong value')
end

cathData = xlsread(cathDataFile,sheet,xlRange);
minSlice = cathData(1,nTimeframe);
maxSlice = cathData(2,nTimeframe);
sliceRange = minSlice:maxSlice; % 35:55
scrSz = get(0, 'Screensize');
featuresAll = struct([]);
featuresAllNorm = struct([]);

%% Reading the data
filename = 'LV Catheter 07.nrrd'; % delete after getting  features
[X, meta] = nrrdread(filename);
sz = sscanf(meta.sizes, '%d');
nDims = sscanf(meta.dimension, '%d');
I = squeeze(X(:,:,:,nTimeframe));

%% Getting the data
% Binarization
for nSlice = sliceRange
    img = GetImage(I, nSlice, ax);
    imshow(img, 'InitialMagnification', 'fit');
    % level = threshTool(img)/255;
    % [level,EM] = graythresh(img);
    [level, EM] = GetThresholdingLevel(img, 'ImprovedOtsu');
    BW = imbinarize(img, level);
    imshow(BW, 'InitialMagnification', 'fit');
    BW = bwareaopen(BW, openArea);
    imshow(BW, 'InitialMagnification', 'fit');
    if isSeparate == 1
        BW = GetSeparatedRegions(BW, method, ax);
    end
    imshow(BW, 'InitialMagnification', 'fit');
    se = strel('ball', 2, 0, 0);
    BW = imopen(BW, se);
    imshow(BW, 'InitialMagnification', 'fit');
    if isVisual == 0
        str1 = sprintf('Binarized image');
        str2 = sprintf('Thresholding level: %.3f (%d out of 255)', level, uint8(level*255));
        str3 = sprintf('Effectiveness metric: %.3f', EM);
        imshow(BW, 'InitialMagnification', 'fit'); AddTitle({str1; str2; str3});
        set(gcf, 'Position', scrSz, 'Color', 'w');
    end
    vars.binarization = {'str1', 'str2', 'str3', 'EM'};
    clear(vars.binarization{:});

    %% Feature extraction
    featuresTemp = struct('Presence', [], ...
                          'Area', [], ... 
                          'ConvexArea', [], ...
                          'Perimeter', [], ...
                          'Eccentricity', [], ...
                          'Solidity', [], ...
                          'Extent', [], ...
                          'EquivDiameter', [], ...
                          'MaxIntensity', [], ...
                          'MeanIntensity', [], ...
                          'Variance', [], ...
                          'StandardDeviation', [], ...
                          'Skewness', [], ...
                          'Kurtosis', [], ...
                          'Contrast', [], ...
                          'Correlation', [], ...
                          'Entropy', [], ...
                          'Energy', [], ...
                          'Homogeneity', [], ...
                          'Centroid', [], ...
                          'Distance', []);

    %% Filling the holes
    if isFill == 1
        BWfill = imfill(BW, 'holes');
    elseif isFill == 0
        BWfill = BW;
    end
    CCfill = bwconncomp(BWfill);
    Lfill = labelmatrix(CCfill);
    numFill = CCfill.NumObjects;
        
    if isVisual == 0
        colorLabelIni = label2rgb(Lfill, colormap(brewermap([],'Accent')), 'k', 'shuffle');
        str1 = sprintf('Region labeling');
        str2 = sprintf('Objects found: %d', numFill);
        imshow(colorLabelIni, 'InitialMagnification', 'fit'); AddTitle({str1; str2});
        set(gcf, 'Position', scrSz, 'Color', 'w');
    end
    vars.fillingHoles = {'str1', 'str2'};
    clear(vars.fillingHoles{:});
    
    %% Main feature analysis
    featuresFill = regionprops(Lfill,  img, {'Area', ... 
                                             'PixelValues', ...
                                             'Eccentricity', ...
                                             'Solidity', ...
                                             'Extent', ...
                                             'MajorAxisLength', ...
                                             'MinorAxisLength', ... 
                                             'EquivDiameter', ...
                                             'MaxIntensity', ...
                                             'MeanIntensity', ...
                                             'Extrema', ...
                                             'BoundingBox',...
                                             'Perimeter', ...
                                             'ConvexArea', ...
                                             'Centroid'});   
    numFill = numel(featuresFill);
    for count = 1:numFill
            featuresFill(count).StandardDeviation = std(double(featuresFill(count).PixelValues));
            featuresFill(count).Variance = var(double(featuresFill(count).PixelValues));
            featuresFill(count).Entropy = entropy(featuresFill(count).PixelValues);
            featuresFill(count).Skewness = skewness(double(featuresFill(count).PixelValues));
            featuresFill(count).Kurtosis = kurtosis(double(featuresFill(count).PixelValues));
    end
    if isVisual == 0
        imshow(BWfill, 'InitialMagnification', 'fit');
        AddTitle('Mean and Standard deviation of regions');
        pmSymbol = char(177);
        hold on
        for count = 1:numFill
            rectangle('Position', featuresFill(count).BoundingBox, 'EdgeColor','c');
            posX = featuresFill(count).Extrema(1,1) - 4;
            posY = featuresFill(count).Extrema(1,2) - 3;
            str = sprintf("%d: %d%s%d", count, round(featuresFill(count).MeanIntensity,0), pmSymbol, round(featuresFill(count).StandardDeviation, 0));
            text(posX, posY, char(str), 'FontSize', 14, 'FontName', 'Times New Roman', 'Color', 'g');
        end
        set(gcf, 'Position', scrSz, 'Color', 'w');
        hold off
    end
    vars.mainFeatureAnalysis = {'count', 'PosX', 'PosY', 'str', 'pmSymbol'};
    clear(vars.mainFeatureAnalysis{:});
    
    %% GLCM analysis
    if ~isempty(featuresFill)
        glcm = cell(1, numFill);
        glcmprops = struct('Contrast', [], ...
                           'Correlation', [], ...
                           'Energy', [], ...
                           'Homogeneity', []);
        for count = 1:numFill
            rect = featuresFill(count).BoundingBox;
            croppedImg = imcrop(img, rect);
            glcm{count} = graycomatrix(croppedImg, 'NumLevels', 255);
            glcmprops(count) = graycoprops(glcm{count}); % all features
%             tempVectorImg = double(croppedImg);
%             tempVectorImg = tempVectorImg(:);
%             featuresFill(count).Skewness = skewness(tempVectorImg); % Compute skewness based on cropped image
%             featuresFill(count).Kurtosis = kurtosis(tempVectorImg); % Compute kurtosis based on cropped image
        end
    else
        glcmprops = struct([]);
    end
    vars.glcmAnalysis = {'i', 'rect', 'str1', 'str2', 'croppedImg', 'colorScheme', ...
                         'posX', 'posY', 'str', 'tempVectorImg'};
    clear(vars.glcmAnalysis{:});
    
    %% Merging all features into one structure 
    numFeatures = numel(featuresFill);
    for count = 1:numFeatures
        for fn = fieldnames(featuresFill)'
           featuresTemp(count).(fn{1}) = featuresFill(count).(fn{1});
        end
    end

    for count = 1:numFeatures
        for fn = fieldnames(glcmprops)'
           featuresTemp(count).(fn{1}) = glcmprops(count).(fn{1});
        end
    end
    
    if isVisual == 1
        hFig = figure;
        imshowpair(BWfill, img, 'montage');
        str0 = sprintf("%d timeframe: %d:%d", nTimeframe,minSlice,maxSlice);
        str1 = sprintf("%d slice", nSlice);
        str2 = 'Enumeration of regions';
        AddTitle({str0, str1, str2});
        hold on
        for count = 1:numFill
            posX = featuresFill(count).Extrema(1,1) - 4;
            posY = featuresFill(count).Extrema(1,2) - 7;
            str = sprintf("%d", count);
            text(posX, posY, char(str), 'FontSize', 18, 'FontName', 'Times New Roman', 'Color', 'g');
%             rectangle('Position', featuresFill(count).BoundingBox, 'EdgeColor','c');
        end
        set(gcf, 'Position', [scrSz(3), 0, scrSz(3), scrSz(4)],...
            'Color', 'w', 'name', str1, 'numbertitle', 'off'); % FOR THE SECOND DISPLAY ONLY
%         set(gcf, 'Position', [scrSz(3)/2, scrSz(2), scrSz(3)/2, scrSz(4)],...
%             'Color', 'w', 'name', str1, 'numbertitle', 'off'); % FOR ONE DISPLAY ONLY    
        hold off
    end
    [featuresTemp(1:numFeatures).Presence] = deal(0);
    
    if isLabel == 1
        cathCoord = ginput(1);   
        [length, height] = size(img);
        xCath = cathCoord(1);
        yCath = cathCoord(2);
        if ~(xCath < 0 || yCath < 0 || xCath >= 2*length || yCath >= height)
            for count = 1:numFill
                featuresTemp(count).Distance = pdist([featuresTemp(count).Centroid; cathCoord], 'euclidean');
            end
            [minDistVal, minDistIndex] = min([featuresTemp.Distance]);
            featuresTemp(minDistIndex).Presence = 1;
        end
    end
    
    fieldsToDel = {'Extrema', 'BoundingBox', 'PixelValues', 'Centroid', 'Distance'};
    featuresTemp = rmfield(featuresTemp,fieldsToDel);
    close(hFig);
    disp(str1);
    
    % Removing NaNs
    fn = fieldnames(featuresTemp);
    for i = 1:numel(fn) 
        featuresTemp = featuresTemp(~isnan([featuresTemp.(fn{i})]));
    end
    featuresAll = [featuresTemp, featuresAll];
    vars.removingNaNs = {'fn', 'i'};
    clear(vars.removingNaNs{:});

    % Normalization
    if isNormalize == 1
        featuresTempNorm = RobustNormalization(featuresTemp, 'quantile', 0);
        featuresTempNorm = featuresTempNorm';
        featuresAllNorm = [featuresTempNorm, featuresAllNorm];
    end
   
    vars.allFeatureAnalysis = {'count', 'PosX', 'PosY', 'str1', 'str2'};
    clear(vars.allFeatureAnalysis{:});
end
disp('Done')
% writetable(struct2table(featuresAll), 'features.xlsx')