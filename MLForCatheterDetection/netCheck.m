%% Initial state
clear all; close all; clc;
set(0, 'DefaultFigureWindowStyle', 'normal');

%% Global variables
nTestTimeframe = 14; %choose from 14 to 17
scrSz = get(0, 'Screensize');
isVisual = 0;

%% Reading the data
tic;
filename = 'LV Catheter 07.nrrd';
[X, meta] = nrrdread(filename); %with meta data
sz = sscanf(meta.sizes, '%d');
I = squeeze(X(:,:,:,nTestTimeframe));
% implay(I);
timerVal = toc;
fprintf('Reading time is %.3f seconds \n', timerVal)
vars.reading = {'filename', 'meta'};
clear(vars.reading{:});


%% Binarization
tic;
nSlice = 56; %1:208
img = I(:,:,nSlice);
% level = threshTool(img)/255;
[level,EM] = graythresh(img);
BW = imbinarize(img, level);
timerVal = toc;
if isVisual == 1
    str1 = sprintf('Binarized image');
    str2 = sprintf('Thresholding level: %.3f (%d out of 255)', level, uint8(level*255));
    str3 = sprintf('Effectiveness metric: %.3f', EM);
    imshow(BW, 'InitialMagnification', 'fit'); addTitle({str1; str2; str3});
    set(gcf, 'Position', scrSz, 'Color', 'w');
end
fprintf('Binarization time is %.4f seconds \n', timerVal)
vars.binarization = {'str1', 'str2', 'str3', 'EM'};
clear(vars.binarization{:});

%% Feature extraction              
featuresReference = struct('Area', [], ...
                     'ConvexArea', [], ...
                     'Perimeter', [], ...
                     'Eccentricity', [], ...
                     'EquivDiameter', [], ...
                     'MaxIntensity', [], ...
                     'MeanIntensity', [], ...
                     'MinIntensity', [], ...
                     'Variance', [], ...
                     'StandardDeviation', [], ...
                     'Contrast', [], ...
                     'Correlation', [], ...
                     'Energy', [], ...
                     'Homogeneity', [], ...
                     'MajorAxisLength', [], ...
                     'MinorAxisLength', [], ...
                     'Orientation', []);

% Extract shape and intensity features
tic;
featuresExtraction = regionprops(BW,  img, {'Area', ...
                                            'ConvexArea', ...
                                            'Perimeter', ...
                                            'Eccentricity', ...
                                            'EquivDiameter', ...
                                            'MaxIntensity', ...
                                            'MeanIntensity', ...
                                            'MinIntensity', ...
                                            'MajorAxisLength', ...
                                            'MinorAxisLength', ...
                                            'Orientation'});
                                        
featuresAdditional = regionprops(BW,  img, {'PixelValues',...
                                            'BoundingBox'});
timerVal = toc;
fprintf('1. Shape features extraction time is %.4f seconds \n', timerVal)
                                  
% Extract STD and variance 
tic;
numExtraction = numel(featuresExtraction);
for count = 1:numExtraction
    featuresExtraction(count).StandardDeviation = std(double(featuresAdditional(count).PixelValues));
    featuresExtraction(count).Variance = var(double(featuresAdditional(count).PixelValues));
end
timerVal = toc;
fprintf('2. STD and variance extraction time is %.4f seconds \n', timerVal)

% Extract texture features
tic;
for count = 1:numExtraction
    rect = featuresAdditional(count).BoundingBox;
    croppedImg = imcrop(img, rect);
    glcm = graycomatrix(croppedImg, 'NumLevels', 255);
    stats = graycoprops(glcm, {'Contrast', 'Correlation', 'Energy', 'Homogeneity'});
    featuresExtraction(count).Contrast = getfield(stats, 'Contrast');
    featuresExtraction(count).Correlation = getfield(stats, 'Correlation');
    featuresExtraction(count).Energy = getfield(stats, 'Energy');
    featuresExtraction(count).Homogeneity = getfield(stats, 'Homogeneity');
end
featuresExtraction = orderfields(featuresExtraction, featuresReference);
timerVal = toc;
fprintf('3. Texture features extraction time is %.4f seconds \n', timerVal)
vars.glcmAnalysis = {'rect', 'croppedImg', 'glcm', 'stats'};
clear(vars.glcmAnalysis{:});

%% Net implementation




