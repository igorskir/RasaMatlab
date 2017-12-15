function [image, count] = DetectMSERblobs(BW)
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

%% Extract MSER Features from an Image (Maximally Stable Extremal Regions)
[featuresMSER,CCmser] = detectMSERFeatures(BW, 'RegionAreaRange', [limits.minArea limits.maxArea]);

%% Eccentricity detection
featuresEccenAna = regionprops(CCmser, {'Area', 'Eccentricity', 'BoundingBox', 'Centroid','Extrema'});
indexEccenAna = find([featuresEccenAna.Eccentricity] >= limits.minEccen & [featuresEccenAna.Eccentricity] <= limits.maxEccen);
featuresEccenDet = featuresMSER(indexEccenAna);
Leccen = labelmatrix(CCmser);
BWeccen = ismember(Leccen, indexEccenAna);
count = featuresEccenDet.Count;
image = BWeccen;

end

