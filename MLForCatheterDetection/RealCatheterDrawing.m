function BWref = RealCatheterDrawing(BW)
sz = size(BW);
[featuresMSER, CCmser] = detectMSERFeatures(BW);
stats = regionprops(CCmser, {'Area', 'Eccentricity', 'BoundingBox', 'Centroid','Extrema'});
numMSER = featuresMSER.Count;
isRefrenceCath = 1; 
a = 9.15; % 16
b = 8.95 ; % 15

% figure;
% imshowpair(BW, 'InitialMagnification', 'Fit');
% hold on;
% plot(featuresMSER,'showPixelList', false,'showEllipses',true)

%% Ellipses drawing
BWroi = zeros(sz(1), sz(2));
BWroi = logical(BWroi);
featuresEllipse = struct('CenterX', [], 'CenterY', [],...
                         'Length', [], 'Width', [],...
                         'Theta',[]);
[Xmesh, Ymesh] = meshgrid(1:sz(1), 1:sz(2));        



% Initial data
featuresEllipse.CenterX = featuresMSER.Location(numMSER,1); % Coordinate X
featuresEllipse.CenterY = featuresMSER.Location(numMSER,2); % Coordinate Y
if isRefrenceCath == 1
    featuresEllipse.Length = b;
    featuresEllipse.Width = a;
else
    featuresEllipse.Length = featuresMSER.Axes(numMSER,1);
    featuresEllipse.Width = featuresMSER.Axes(numMSER,2);
end
featuresEllipse.Theta = rad2deg(featuresMSER.Orientation(numMSER)); % degree you want to rotate

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
centerOrigX = stats.Centroid(1);
centerOrigY = stats.Centroid(2);
widthNew = size(rotImg,1);
lenghtNew = size(rotImg,2);
PosX = centerOrigX - lenghtNew/2; % 2 version
PosY = centerOrigY - widthNew/2; % 2 version
offset = 0;
PosX = round(PosX,0) - offset;
PosY = round(PosY,0) - offset;
BWroi((1:size(rotImg,1)) + PosY,(1:size(rotImg,2)) + PosX,:) = rotImg;


BWref = BWroi;

end

