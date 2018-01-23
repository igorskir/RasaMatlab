%% Initial state of the programm
clear all; close all; clc;
set(0, 'DefaultFigureWindowStyle', 'normal');
currentFolder = pwd;
addpath(genpath(pwd));
% opengl('save', 'software')
% opengl hardware

%% Limits
limits = struct('maxBlobsArea', 0,... % 5
                'minArea', 65, 'maxArea', 195,... % 5 200    % USED
                'minEccen', 0.342345237561179, 'maxEccen', 0.899968089081084,...   % 0 0.96    % USED
                'minStd', 3, 'maxStd', 53,...     % 3 53
                'minMean', 60, 'maxMean', 132,... % 74 132
                'minContrast', 45, 'maxContrast', 574,... % 77 574
                'minCorrelation', 0.55, 'maxCorrelation', 0.97,... % 0.55 0.97
                'minEnergy', 0.002, 'maxEnergy', 0.112,... % 0.002 0.112
                'minHomogeneity', 0.07, 'maxHomogeneity', 0.27,... % 0.07 0.27
                'minVolume', 2600, 'maxVolume', 10000); % 300 10000 

%% Global variables
isVisual = 1;
nTimeframe = 9; %9
nSlice = 80; %60
scrSz = get(0, 'Screensize');

% % for i = 67:109 % FOR MEASUREMENT OF OBJECT NUMBER 
% % nSlice = i;

%% Reading the data
tic;
filename = 'LV Catheter 07.nrrd'; % delete after getting  features
[X, meta] = nrrdread(filename);
Y = double(X);
sz = sscanf(meta.sizes, '%d');
nDims = sscanf(meta.dimension, '%d');
toc;

%% Binarization
tic;
I = squeeze(X(:,:,:,nTimeframe)); % short-axis view
img = I(:,:,nSlice); % short-axis view
% img = I(:,nSlice,:);   % long-axis view
% img = squeeze(img);    % long-axis view
% img = rot90(img);      % long-axis view
% level = threshTool(img)/255;
[level,EM] = graythresh(img);
% level = 0.2057;
% EM = 1;
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

%% Filling holes
BWfill = imfill(BW, 'holes');
CCfill = bwconncomp(BWfill);
Lfill = labelmatrix(CCfill);
numFill = CCfill.NumObjects;

if isVisual == 1
    colorLabelFill = label2rgb(Lfill, 'parula', 'k', 'shuffle');
    str1 = sprintf('Region labeling');
    str2 = sprintf('Objects found: %d', numFill);
    imshow(colorLabelFill, 'InitialMagnification', 'fit'); addTitle({str1; str2});
    set(gcf, 'Position', scrSz, 'Color', 'w');
    vars.fillingHoles = {'str1', 'str2'};
    clear(vars.fillingHoles{:});
end

%% Enumeration of the regions
if isVisual == 1
    str1 = sprintf('Region enumeration');
    str2 = sprintf('Objects found: %d', numFill);
    vislabels(Lfill); addTitle({str1; str2});
    set(gcf, 'Position', scrSz, 'Color', 'w');
    vars.enumeration = {'str1', 'str2'};
    clear(vars.enumeration{:});
end

%% Extract MSER Features from an Image (Maximally Stable Extremal Regions)
[featuresMSER,CCmser] = detectMSERFeatures(BWfill, 'RegionAreaRange', [limits.minArea limits.maxArea]);
numMSER = featuresMSER.Count;
if isVisual == 1
    figure('Position', [scrSz(1), scrSz(2), scrSz(3)/2, scrSz(4)]); 
    imshow(BWfill, 'InitialMagnification', 'Fit');
    str1 = sprintf('Extract MSER features');
    str2 = sprintf('Objects found: %d', numMSER);
    addTitle({str1, str2});
    hold on;
    plot(featuresMSER);
    figure('Position', [scrSz(3)/2, scrSz(2), scrSz(3)/2, scrSz(4)]); 
    imshow(BWfill, 'InitialMagnification', 'Fit');
    addTitle({str1, str2});
    hold on;
    plot(featuresMSER, 'showPixelList', true, 'showEllipses', false);
    vars.extractMSER = {'str1', 'str2'};
    clear(vars.extractMSER{:});
end

%% Eccentricity detection
featuresEccenAna = regionprops(CCmser, {'Area', 'Eccentricity', 'BoundingBox', 'Centroid','Extrema'});
indexEccenAna = find([featuresEccenAna.Eccentricity] >= limits.minEccen & [featuresEccenAna.Eccentricity] <= limits.maxEccen);
featuresEccenDet = featuresMSER(indexEccenAna);
Leccen = labelmatrix(CCmser);
BWeccen = ismember(Leccen, indexEccenAna);
numEccen = featuresEccenDet.Count;

if isVisual == 1
    figure('Position', [scrSz(1), scrSz(2), scrSz(3)/2, scrSz(4)], 'Color', 'w'); 
    imshow(BWfill, 'InitialMagnification', 'fit');
    str1 = sprintf('Eccentricity of the regions');
    str2 = sprintf('Objects found: %d', numMSER);
    addTitle({str1, str2});
    hold on
    for count = 1:numMSER
        rectangle('Position', featuresEccenAna(count).BoundingBox, 'EdgeColor','c');
        posX = featuresEccenAna(count).Extrema(1,1) - 8;
        posY = featuresEccenAna(count).Extrema(1,2) - 5;
        str3 = sprintf("%d: %0.2f", count, round(featuresEccenAna(count).Eccentricity, 2));
        text(posX, posY, char(str3), 'FontSize', 32, 'FontName', 'Times New Roman', 'Color', 'g');
    end

    figure('Position', [scrSz(3)/2, scrSz(2), scrSz(3)/2, scrSz(4)], 'Color', 'w'); 
    imshow(BWfill, 'InitialMagnification', 'Fit');
    str4 = sprintf('Eccentricity detection');
    str5 = sprintf('Objects found: %d', numEccen);
    addTitle({str4, str5});
    hold on;
    plot(featuresEccenDet,'showPixelList',true,'showEllipses',false);
    vars.detectEccentricity = {'str1', 'str2', 'str3', 'str4', 'str5'};
    clear(vars.detectEccentricity{:});
end
toc;
%% Ellipses drawing
BWroi = zeros(sz(1), sz(2));
BWroi = logical(BWroi);
isVisual = 0;
featuresEllipse = struct('CenterX', [], 'CenterY', [],...
                         'Length', [], 'Width', [],...
                         'Theta',[]);
[Xmesh, Ymesh] = meshgrid(1:sz(1), 1:sz(2));        
for count = 1:numEccen % 1
    % Initial data
    featuresEllipse.CenterX = featuresEccenDet.Location(count,1); % Coordinate X
    featuresEllipse.CenterY = featuresEccenDet.Location(count,2); % Coordinate Y
    featuresEllipse.Length = featuresEccenDet.Axes(count,1);
    featuresEllipse.Width = featuresEccenDet.Axes(count,2);
    featuresEllipse.Theta = rad2deg(featuresEccenDet.Orientation(count)); % degree you want to rotate
    
    % Binary Mask which multiplies to an image
    ellipse = ((Xmesh - featuresEllipse.CenterX)/(0.5*featuresEllipse.Length)).^2 + ...
             +((Ymesh - featuresEllipse.CenterY)/(0.5*featuresEllipse.Width)).^2 <= 1;
    if isVisual == 1
        imshow(ellipse);
    end
    
    % Bounding box for an ellipse
    CCreg = bwconncomp(ellipse);
    Lreg = labelmatrix(CCreg);
    featuresROI = regionprops(Lreg, {'BoundingBox', 'Centroid'});
    numROI = numel(featuresROI);
    if isVisual == 1
        imshow(ellipse, 'InitialMagnification', 'Fit');
        hold on
        for count = 1:numROI
            rectangle('Position', featuresROI(count).BoundingBox, 'EdgeColor', 'c');
        end
        hold off
    end
    
    % Crop
    rect = featuresROI(numROI).BoundingBox;
    rect(3) = rect(3) - 1;
    rect(4) = rect(4) - 1;
    cropImg = imcrop(ellipse, rect);
    if isVisual == 1
        imtool(cropImg);
    end
    
    % Rotation
    rotImg = imrotate(cropImg, featuresEllipse.Theta, 'bilinear');
    if isVisual == 1
        imtool(rotImg);
    end
    
    % Superimpose
    centerOrigX = featuresEccenAna(indexEccenAna(count)).Centroid(1);
    centerOrigY = featuresEccenAna(indexEccenAna(count)).Centroid(2);
    widthNew = size(rotImg,1);
    lenghtNew = size(rotImg,2);
    PosX = centerOrigX - lenghtNew/2; % 2 version
    PosY = centerOrigY - widthNew/2; % 2 version
    offset = 1;
    PosX = round(PosX,0) - offset;
    PosY = round(PosY,0) - offset;
    BWroi((1:size(rotImg,1)) + PosY,(1:size(rotImg,2)) + PosX,:) = rotImg;
end
%%
% isVisual = ~isVisual;
if isVisual == 0
    coloredROI = label2rgb(BWroi, 'parula', 'k', 'shuffle');
    fusedImg = imfuse(BWfill,coloredROI, 'blend');
    imshow(fusedImg, 'InitialMagnification', 'Fit');
    hold on;
    plot(featuresEccenDet);
end

% end % FOR NUMBER OF OBJECTS MEASUREMENT