%%
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
%     PosX = round(featuresEccenAna(indexEccenAna(count)).BoundingBox(1), 0); % 1 version
%     PosY = round(featuresEccenAna(indexEccenAna(count)).BoundingBox(2), 0); % 1 version
    PosX = centerOrigX - lenghtNew/2; % 2 version
    PosY = centerOrigY - widthNew/2; % 2 version
    offset = 1;
    PosX = round(PosX,0) - offset;
    PosY = round(PosY,0) - offset;
    BWroi((1:size(rotImg,1)) + PosY,(1:size(rotImg,2)) + PosX,:) = rotImg;
end

if isVisual == 0
    coloredROI = label2rgb(BWroi, 'parula', 'k', 'shuffle');
    fusedImg = imfuse(BWfill,coloredROI, 'blend');
    imshow(fusedImg, 'InitialMagnification', 'Fit');
    hold on;
    plot(featuresEccenDet);
end

