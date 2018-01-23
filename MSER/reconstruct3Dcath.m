%% Initialization
tic;
isEllipse = 1; % draw ellipses or not   %%UNCOMMENT AFTER TIME TESTING
% % % Main data
% % filename = 'LV Catheter 07.nrrd';
% % [X, meta] = nrrdread(filename);
% % sz = sscanf(meta.sizes, '%d');
% % nDims = sscanf(meta.dimension, '%d');
% % nTimeframe = 9;
% % I = squeeze(X(:,:,:,nTimeframe));
% % BWr = zeros(sz(1), sz(2), sz(3));

% Ellipse data
[Xmesh, Ymesh] = meshgrid(1:sz(1), 1:sz(2));
BWroi = zeros(sz(1), sz(2));
BWroi = logical(BWroi);
featuresEllipse = struct('CenterX', [], 'CenterY', [],...
                         'Length', [], 'Width', [],...
                         'Theta',[]);
                    
for nSlice = 1:sz(3) 
    img = I(:,:,nSlice);
    % Binarization
%     level = graythresh(img); % MAIN VERSION
    level = graythresh(img) + 0.1*graythresh(img); % Try to increase accuracy by incresaing threshold
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
    featuresEccenAna = regionprops(CCmser, {'Area', 'Eccentricity', 'BoundingBox', 'Centroid','Extrema'});
    indexEccenAna = find([featuresEccenAna.Eccentricity] >= limits.minEccen & [featuresEccenAna.Eccentricity] <= limits.maxEccen);
    featuresEccenDet = featuresMSER(indexEccenAna);
    Leccen = labelmatrix(CCmser);
    BWeccen = ismember(Leccen, indexEccenAna);
    numEccen = featuresEccenDet.Count;
    % Draw ellipses
    if isEllipse == 1
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
    else
        BWroi = BWeccen;   
    end
    BWr(:,:,nSlice) = BWroi;
    BWroi = zeros(sz(1), sz(2));
    BWroi = logical(BWroi);
end

% Volume detection
if ~isempty(BWr)
    CCvol = bwconncomp(BWr);
    Lvol = labelmatrix(CCvol);
    featuresVol = regionprops(CCvol,'Centroid','PixelIdxList');
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
toc;
% showvol(BWs, 0.1)
% implay(BWr)