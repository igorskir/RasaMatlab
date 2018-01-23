%% Initialization
tic;
clear all; close all; clc;
set(0, 'DefaultFigureWindowStyle', 'normal');
currentFolder = pwd;
addpath(genpath(pwd));

% Main data
filename = 'LV Catheter 07.nrrd';
[X, meta] = nrrdread(filename);
sz = sscanf(meta.sizes, '%d');
nDims = sscanf(meta.dimension, '%d');

limits = struct('maxBlobsArea', 0,... % 5
                'minArea', 1, 'maxArea', 200,... % 5 200
                'minEccen', 0.0, 'maxEccen', 0.96,...
                'minStd', 3, 'maxStd', 53,...     % 3 53
                'minMean', 60, 'maxMean', 132,... % 74 132
                'minContrast', 45, 'maxContrast', 574,... % 77 574
                'minCorrelation', 0.55, 'maxCorrelation', 0.97,... % 0.55 0.97
                'minEnergy', 0.002, 'maxEnergy', 0.112,... % 0.002 0.112
                'minHomogeneity', 0.07, 'maxHomogeneity', 0.27,... % 0.07 0.27
                'minVolume', 1900, 'maxVolume', 10000); % 300 10000

% Ellipse data
[Xmesh, Ymesh] = meshgrid(1:sz(1), 1:sz(2));
BWroi = zeros(sz(1), sz(2));
BWroi = logical(BWroi);
featuresEllipse = struct('CenterX', [], 'CenterY', [],...
                         'Length', [], 'Width', [],...
                         'Theta',[]);
BWfull = zeros(sz(1), sz(2), sz(3), sz(4));
objectInfo = zeros(1,5);
for nTimeframe = 1:sz(4)
    BWr = zeros(sz(1), sz(2), sz(3));
    I = squeeze(X(:,:,:,nTimeframe));
    for nSlice = 1:sz(3) 
        img = I(:,:,nSlice);
        % Binarization
%         level = graythresh(img);
        level = 0.1557;
        BW = imbinarize(img, level);
        % Filling holes
        BWfill = imfill(BW, 'holes');
        CCfill = bwconncomp(BWfill);
        Lfill = labelmatrix(CCfill);
        numFill = CCfill.NumObjects; 
        % Extract MSER Features from an Image (Maximally Stable Extremal Regions)
        [featuresMSER,CCmser] = detectMSERFeatures(BWfill, 'RegionAreaRange', [limits.minArea limits.maxArea]);
        numMSER = featuresMSER.Count;
        % Eccentricity detection
        featuresEccenAna = regionprops(CCmser, {'Area', 'Eccentricity', 'BoundingBox', 'Centroid', 'Extrema'});
        indexEccenAna = find([featuresEccenAna.Eccentricity] >= limits.minEccen & [featuresEccenAna.Eccentricity] <= limits.maxEccen);
        featuresEccenDet = featuresMSER(indexEccenAna);
        numEccen = featuresEccenDet.Count;
        % Draw ellipses
        for count = 1:numEccen   
        % Initial data
            featuresEllipse.CenterX = featuresEccenDet.Location(count,1); % Coordinate X
            featuresEllipse.CenterY = featuresEccenDet.Location(count,2); % Coordinate Y
            featuresEllipse.Length = featuresEccenDet.Axes(count,1);
            featuresEllipse.Width = featuresEccenDet.Axes(count,2);
            featuresEllipse.Theta = rad2deg(featuresEccenDet.Orientation(count)); % degree you want to rotate

            % Binary Mask which multiplies to an image
            ellipse = ((Xmesh - featuresEllipse.CenterX)/(0.5*featuresEllipse.Length)).^2 + ...
                     +((Ymesh - featuresEllipse.CenterY)/(0.5*featuresEllipse.Width)).^2 <= 1;

            % Bounding box for an ellipse
            CCreg = bwconncomp(ellipse);
            Lreg = labelmatrix(CCreg);
            featuresROI = regionprops(Lreg, {'BoundingBox', 'Centroid'});
            numROI = numel(featuresROI);

            % Crop
            rect = featuresROI(numROI).BoundingBox;
            rect(3) = rect(3) - 1;
            rect(4) = rect(4) - 1;
            cropImg = imcrop(ellipse, rect);

            % Rotation
            rotImg = imrotate(cropImg, featuresEllipse.Theta, 'bilinear');

            % Superimpose
            centerOrigX = featuresEccenAna(indexEccenAna(count)).Centroid(1);
            centerOrigY = featuresEccenAna(indexEccenAna(count)).Centroid(2);
            widthNew = size(rotImg,1);
            lenghtNew = size(rotImg,2);
        %     PosX = round(featuresEccenAna(indexEccenAna(count)).BoundingBox(1), 0); % 1 version
        %     PosY = round(featuresEccenAna(indexEccenAna(count)).BoundingBox(2), 0); % 1 version
            PosX = centerOrigX - lenghtNew/2; % 2 version
            PosY = centerOrigY - widthNew/2; % 2 version
            offset = 0;
            PosX = round(PosX,0) - offset;
            PosY = round(PosY,0) - offset;
            BWroi((1:size(rotImg,1)) + PosY,(1:size(rotImg,2)) + PosX,:) = rotImg;
        end
        BWr(:,:,nSlice) = BWroi;
        CCr = bwconncomp(BWr);
        BWroi = zeros(sz(1), sz(2));
        BWroi = logical(BWroi);
    end

    % Volume detection
    if ~isempty(BWr)
        CCvol = bwconncomp(BWr);
        Lvol = labelmatrix(CCvol);
        featuresVol = regionprops(CCvol, 'Centroid', 'PixelIdxList');
        numVol = numel(featuresVol);
        for idx = 1:numVol
            featuresVol(idx).Volume = numel(featuresVol(idx).PixelIdxList);
        end
        index = find([featuresVol.Volume] > limits.minVolume & [featuresVol.Volume] < limits.maxVolume);
        BWvol = ismember(Lvol, index);
        BWvol = double(BWvol);
    else
        BWvol = zeros();
    end
    BWs = smooth3(BWvol, 'box', [3,3,5]);
    CCs = bwconncomp(BWs);
    CCvol = bwconncomp(BWvol);
    objectInfo(nTimeframe,1) = nTimeframe;
    objectInfo(nTimeframe,2) = CCr.NumObjects;
    objectInfo(nTimeframe,3) = CCvol.NumObjects;
    objectInfo(nTimeframe,4) = CCs.NumObjects;
    for i = 1:CCr.NumObjects
        [temp , ~] = size(CCr.PixelIdxList{1,i});
        objectVolume(1,i) = temp;
    end
    objectInfo(nTimeframe,5) = min(objectVolume);
    objectInfo(nTimeframe,6) = max(objectVolume);
    BWfull(:,:,:,nTimeframe) = BWs;
end
toc;