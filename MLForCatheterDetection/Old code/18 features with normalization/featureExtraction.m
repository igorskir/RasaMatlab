%% Initial state of the programm
clear all; close all; clc;
set(0, 'DefaultFigureWindowStyle', 'normal');
currentFolder = pwd;
addpath(genpath(pwd));

%% Global variables
isNormalize = 1;
isVisual = 0;
isFill = 0;
nTimeframe = 1; %9
%XLS reading
cathDataFile = 'D:\Clouds\Google Drive\RASA\Matlab\12. Catheter ML\Presentation\LV Catheter 07.xlsx';
sheet = 1;
xlRange = 'C7:S8';
cathData = xlsread(cathDataFile,sheet,xlRange);
minSlice = cathData(1,nTimeframe);
maxSlice = cathData(2,nTimeframe);
sliceRange = minSlice:maxSlice; % 35:55
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
    if isVisual == 1
        str1 = sprintf('Binarized image');
        str2 = sprintf('Thresholding level: %.3f (%d out of 255)', level, uint8(level*255));
        str3 = sprintf('Effectiveness metric: %.3f', EM);
        imshow(BW, 'InitialMagnification', 'fit'); addTitle({str1; str2; str3});
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
                          'MinIntensity', [], ...
                          'Variance', [], ...
                          'StandardDeviation', [], ...
                          'Contrast', [], ...
                          'Correlation', [], ...
                          'Energy', [], ...
                          'Homogeneity', []);

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
    %% Main feature analysis
    tic;
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
                                             'MinIntensity', ...
                                             'Extrema', ...
                                             'BoundingBox',...
                                             'Perimeter', ...
                                             'ConvexArea'});   
    numFill = numel(featuresFill);
    for count = 1:numFill
            featuresFill(count).StandardDeviation = std(double(featuresFill(count).PixelValues));
            featuresFill(count).Variance = var(double(featuresFill(count).PixelValues));
    end
    if isVisual == 1
        imshow(BWfill, 'InitialMagnification', 'fit');
        addTitle('Mean and Standard deviation of regions');
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
        glcmprops = struct('Contrast', [], 'Correlation', [], 'Energy', [], 'Homogeneity', []);
        for count = 1:numFill
            rect = featuresFill(count).BoundingBox;
            croppedImg = imcrop(img, rect);
            glcm{count} = graycomatrix(croppedImg, 'NumLevels', 255);
            glcmprops(count) = graycoprops(glcm{count});
        end
    else
        glcmprops = struct([]);
    end
    vars.glcmAnalysis = {'i', 'rect', 'str1', 'str2', 'croppedImg', 'colorScheme', 'posX', 'posY', 'str'};
    clear(vars.glcmAnalysis{:});
    %% All features
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
    if isVisual == 0
        hFig = figure;
%         hFrame = get(handle(gcf),'JavaFrame'); % MAXIMIZE
%         imshow(BWfill, 'InitialMagnification', 'fit');
        imshowpair(BWfill, img, 'montage');
        str0 = sprintf("%d timeframe: %d:%d", nTimeframe,minSlice,maxSlice);
        str1 = sprintf("%d slice", nSlice);
        str2 = 'Enumeration of regions';
        addTitle({str0, str1, str2});
        hold on
        for count = 1:numFill
            posX = featuresFill(count).Extrema(1,1) - 4;
            posY = featuresFill(count).Extrema(1,2) - 7;
            str = sprintf("%d", count);
            text(posX, posY, char(str), 'FontSize', 18, 'FontName', 'Times New Roman', 'Color', 'g');
        end
%         hFrame.setMaximized(1);  % MAXIMIZE
%         set(gcf, 'Position', scrSz/2, 'Color', 'w', 'name', str1, 'numbertitle', 'off');
% %         set(gcf, 'Position', [scrSz(3), 0, scrSz(3), scrSz(4)],...
% %             'Color', 'w', 'name', str1, 'numbertitle', 'off'); % FOR THE SECOND DISPLAY ONLY
        set(gcf, 'Position', [scrSz(3)/2, scrSz(2), scrSz(3)/2, scrSz(4)],...
            'Color', 'w', 'name', str1, 'numbertitle', 'off'); % FOR ONE DISPLAY ONLY    
%         set(0,'DefaultFigureWindowStyle','docked');
%         set(hFrame,'Maximized',1); % MAXIMIZE
        hold off
    end
    [featuresTemp(1:numFeatures).Presence] = deal(0);
    fieldsToDel = {'Extrema', 'BoundingBox', 'PixelValues'};
    featuresTemp = rmfield(featuresTemp,fieldsToDel);
    
    if isNormalize == 1
        featuresTemp = robustNormalization(featuresTemp, 'std', 0);
        featuresTemp = featuresTemp';
    end
    
    close(hFig);
    disp(str1);
    featuresAll = [featuresTemp, featuresAll];
    vars.allFeatureAnalysis = {'count', 'PosX', 'PosY', 'str1', 'str2'};
    clear(vars.allFeatureAnalysis{:});
end % comment if you need to check 1 slice

% writetable(struct2table(featuresAll), 'features.xlsx')