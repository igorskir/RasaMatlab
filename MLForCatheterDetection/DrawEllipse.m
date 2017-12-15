%% Extract MSER Features from an Image (Maximally Stable Extremal Regions)
img = imread('C:\Users\Viacheslav Danilov\Desktop\img.png');
sz = size(img);
level = graythresh(img);
BW = imbinarize(img, level);
[featuresMSER, CCmser] = detectMSERFeatures(BW);
stats = regionprops(CCmser, {'Area', 'Eccentricity', 'BoundingBox', 'Centroid','Extrema'});
numMSER = featuresMSER.Count;
a = 9.15; % 16
b = 8.95 ; % 15
isVisual = 1;
isRefrenceCath = 1;
scrSz = get(0, 'Screensize');

if isVisual == 1
    imshow(BW, 'InitialMagnification', 'Fit');
    hold on;
    plot(featuresMSER,'showPixelList', false,'showEllipses',true)
end

%% Ellipses drawing
BWroi = zeros(sz(1), sz(2));
BWroi = logical(BWroi);
featuresEllipse = struct('CenterX', [], 'CenterY', [],...
                         'Length', [], 'Width', [],...
                         'Theta',[]);
[Xmesh, Ymesh] = meshgrid(1:sz(1), 1:sz(2));        
for count = 1:numMSER
    % Initial data
    featuresEllipse.CenterX = featuresMSER.Location(count,1); % Coordinate X
    featuresEllipse.CenterY = featuresMSER.Location(count,2); % Coordinate Y
    featuresEllipse.Length = featuresMSER.Axes(count,1);
    featuresEllipse.Width = featuresMSER.Axes(count,2);
    featuresEllipse.Theta = rad2deg(featuresMSER.Orientation(count)); % degree you want to rotate
    
    % Binary Mask which multiplies to an image
    ellipse = ((Xmesh - featuresEllipse.CenterX)/(0.5*featuresEllipse.Length)).^2 + ...
             +((Ymesh - featuresEllipse.CenterY)/(0.5*featuresEllipse.Width)).^2 <= 1;
    if isVisual == 0
        imshow(ellipse);
    end
    
    % Bounding box for an ellipse
    CCreg = bwconncomp(ellipse);
    Lreg = labelmatrix(CCreg);
    featuresROI = regionprops(Lreg, {'BoundingBox', 'Centroid'});
    numROI = numel(featuresROI);
    if isVisual == 0
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
    if isVisual == 0
        imtool(cropImg);
    end
    
    % Rotation
    rotImg = imrotate(cropImg, featuresEllipse.Theta, 'bilinear');
    if isVisual == 0
        imtool(rotImg);
    end
    
    % Superimpose
    centerOrigX = stats.Centroid(1);
    centerOrigY = stats.Centroid(2);
    if isRefrenceCath == 1 
        widthNew = 2*b;
        lenghtNew = 2*a;
    else
        widthNew = size(rotImg,1);
        lenghtNew = size(rotImg,2);
    end
    PosX = centerOrigX - lenghtNew/2; % 2 version
    PosY = centerOrigY - widthNew/2; % 2 version
    offset = 0;
    PosX = round(PosX,0) - offset;
    PosY = round(PosY,0) - offset;
    BWroi((1:size(rotImg,1)) + PosY,(1:size(rotImg,2)) + PosX,:) = rotImg;
end

%% Blending

idxDice = dice(BW, BWroi);
idxJaccard = jaccard(BW, BWroi);
idxBFscore = bfscore(BW, BWroi);

C = imfuse(BW,BWroi, 'ColorChannels', 'red-cyan');
imshow(C, 'InitialMagnification', 'Fit');
str1 = sprintf("Blended overlay image");
str3 = sprintf("Dice index: %.2f, Jaccard index: %.2f, Contour matching score: %.2f", ...
                idxDice, idxJaccard, idxBFscore);
str2 = sprintf("(red for BW, cyan for BWroi, white for similar areas between the two images)");
AddTitle({str1; str2; str3});
set(gcf, 'Position', [scrSz(3), 0, scrSz(3), scrSz(4)],...
    'Color', 'w', 'name', str1, 'numbertitle', 'off'); 








