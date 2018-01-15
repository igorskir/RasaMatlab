%% Initial state
clear all;  close all; clc;
set(0, 'DefaultFigureWindowStyle', 'normal');
addpath(genpath(pwd));

% Settings
net = 'ffb';                    % ffb, cfb, r
filename = '7 timeframe.mat';   % '7 timeframe.mat' or '17 timeframe.mat'
ax = 'short';                   % short, long1, long2
isFill = 0;                     % filling the holes in the objects (1) or not(0)
isSeparate = 0;                 % perform separation of catheter and pigtail (1) or not (0)
isNormalize = 1;                % normolized the data within the column (1) or not (0)
isOverlaySegmentation = 1;      % superimpose catheter segmentation (1) or not (0)
isOverlayBbox = 1;              % superimpose catheter bounding box (1) or not (0)
isSaveImages = 1;               % saving images into the folder (1, NOTE: isVisual = 1) or not (0)        
isVisual = 1;                   % visualization on (1) and off (0)
openArea = 15;
scrSz = get(0, 'Screensize');

% Load a timeframe
cd('MAT files\Timeframes');
load(filename);
cd ..\..

% Load net object
cd('Net tests\20 features');
switch net
    case 'ffb'
        load('feed-forward 20_features [40 20 10].mat', 'net'); 
    case 'cfb'
        load('cascade 20_features [20 10 5].mat', 'net'); 
    case 'r'
        load('recurrent 20_features [20 10 5].mat', 'net'); 
    otherwise
        disp('Choose a network properly')
end
cd ..\..
% processingTime = zeros(size(I,3), 9);
for nSlice = 1:size(I,3)
    tic;
    img = GetImage(I, nSlice, ax);
    
    %% Thresholding and mathematical morphology
    [level, ~] = GetThresholdingLevel(img, 'ImprovedOtsu');
    BW = imbinarize(img, level);
    processingTime(nSlice,1) = toc;    
    tic;
    BW = bwareaopen(BW, openArea);
    if isSeparate == 1
        BW = GetSeparatedRegions(BW, method, ax);
    end
    se = strel('ball', 3, 0, 0);
    BW = imopen(BW, se);
    processingTime(nSlice,2) = toc;
    %% Filling the holes
    tic;
    if isFill == 1
        BW = imfill(BW, 'holes');
    end
    
    %% Get label matrix
    CC = bwconncomp(BW);
    L = labelmatrix(CC);
    numObjects = CC.NumObjects;
    
    %% Check the number of binarized objects
    if numObjects < 2
        if isVisual == 1
            imshow(img, 'InitialMagnification', 'fit');
            figureName = strcat(num2str(nSlice), ' slice');
            AddTitle(figureName);
            set(gcf, 'Position', scrSz, 'Color', 'w');
        end
        if isSaveImages == 1          
            cd('Detection results');
            mkdir(filename)
            cd(filename)
            saveas(gcf,num2str(nSlice),'tiff');
            cd ..\..
        end
        continue;
    end
    processingTime(nSlice,3) = toc;
    %% Main feature analysis
    tic;
    statsTemp = struct('Presence', [], ...
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
                       'Centroid', []);
    
    statsMain = regionprops(L,  img, {'Area', ... 
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
    for count = 1:numObjects
            statsMain(count).StandardDeviation = std(double(statsMain(count).PixelValues));
            statsMain(count).Variance = var(double(statsMain(count).PixelValues));
            statsMain(count).Entropy = entropy(statsMain(count).PixelValues);
            statsMain(count).Skewness = skewness(double(statsMain(count).PixelValues));
            statsMain(count).Kurtosis = kurtosis(double(statsMain(count).PixelValues));
    end
    processingTime(nSlice,4) = toc;
    
    %% GLCM analysis
    tic;
    glcm = cell(1, numObjects);
    statsGLCM = struct('Contrast', [], ...
                       'Correlation', [], ...
                       'Energy', [], ...
                       'Homogeneity', []);
    for count = 1:numObjects
        rect = statsMain(count).BoundingBox;
        croppedImg = imcrop(img, rect);
        glcm{count} = graycomatrix(croppedImg, 'NumLevels', 255);
        statsGLCM(count) = graycoprops(glcm{count}); % all features
    end
    processingTime(nSlice,5) = toc;
    
    %% Merging all features into one structure
    tic;
    for count = 1:numObjects
        for fn = fieldnames(statsMain)'
           statsTemp(count).(fn{1}) = statsMain(count).(fn{1});
        end
    end
    
    for count = 1:numObjects
        for fn = fieldnames(statsGLCM)'
           statsTemp(count).(fn{1}) = statsGLCM(count).(fn{1});
        end
    end
    [statsTemp(1:numObjects).Presence] = deal(0);
    fieldsToDel = {'Extrema', 'BoundingBox', 'PixelValues', 'Centroid'};
    statsTest = rmfield(statsTemp,fieldsToDel);
    processingTime(nSlice,6) = toc;
    
    %% Removing NaNs
    fn = fieldnames(statsTest);
    for i = 2:numel(fn) 
        statsTest = statsTest(~isnan([statsTest.(fn{i})]));
    end
       
    %% Normalization
    tic;
    if isNormalize == 1
        statsTestNorm = RobustNormalization(statsTest, 'quantile', 0);
    end
    processingTime(nSlice,7) = toc;
    
    %% Net application
    tic;
    statsTestNorm = struct2mat(statsTestNorm)';
    for i = 1:numObjects
        statsTestNorm(1,i) = round(net(statsTestNorm(2:end,i)));
    end
    statsTestNorm = statsTestNorm';
    processingTime(nSlice,8) = toc;
    
    %% Visualization
    % Bounding box
    if isVisual == 1
        tic;
        [cathPresence, cathLabelVal] = max(statsTestNorm(:,1));
        hold on
        if cathPresence == 1
            if isOverlaySegmentation == 1
                mask = zeros(size(I,1),size(I,2));
                mask(L == cathLabelVal) = 1;
                imshow(img, 'InitialMagnification', 'fit');
                alphamask(mask, [0 1 0], 0.25);
            end
            if isOverlayBbox == 1
                for count = 1:numObjects
                    if statsTestNorm(count,1) == 1
                        rectangle('Position', statsTemp(count).BoundingBox, 'EdgeColor','y');
            %             posX = statsTemp(count).Extrema(1,1) - 4;
            %             posY = statsTemp(count).Extrema(1,2) - 3;
            %             str = sprintf("Catheter");
            %             text(posX, posY, char(str), 'FontSize', 14, 'FontName', 'Times New Roman', 'Color', 'y');          
                    end
                end
            end
        else
            imshow(img, 'InitialMagnification', 'fit');
        end
        figureName = strcat(num2str(nSlice), ' slice');
        AddTitle(figureName);
        set(gcf, 'Position', scrSz, 'Color', 'w');
        hold off
        processingTime(nSlice, 9) = toc;
    end
    % Saving images 
    if isSaveImages == 1
        cd('Detection results');
        mkdir(filename)
        cd(filename)
        saveas(gcf,num2str(nSlice),'tiff');
        cd ..\..
    end
    fprintf('%.d slice processed \n', nSlice);
end