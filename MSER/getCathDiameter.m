%% Automatically take measurements for diameter
path = 'D:\Clouds\Google Drive\RASA\Matlab\11. MSER catheter detection\Edited (110% threshold)\';
scrSz = get(0, 'Screensize');
% cathRange = 53:108; % 100%
cathRange = 54:104;
isVisual = 0;
diameter = zeros(cathRange(end)-cathRange(1),2);
str = 1;
for numSlice = cathRange
    BW = imread([path, num2str(numSlice),'.png']);
    [featuresMSER,CCmser] = detectMSERFeatures(BW);
    numMSER = featuresMSER.Count;
    if isVisual == 1
        figure('Position', [scrSz(1), scrSz(2), scrSz(3)/2, scrSz(4)]); 
        imshow(BW, 'InitialMagnification', 'Fit');
        str1 = sprintf('Extract MSER features');
        str2 = sprintf('Objects found: %d', numMSER);
        addTitle({str1, str2});
        hold on;
        plot(featuresMSER);
        figure('Position', [scrSz(3)/2, scrSz(2), scrSz(3)/2, scrSz(4)]); 
        imshow(BW, 'InitialMagnification', 'Fit');
        addTitle({str1, str2});
        hold on;
        plot(featuresMSER, 'showPixelList', true, 'showEllipses', false);
        vars.extractMSER = {'str1', 'str2'};
        clear(vars.extractMSER{:});
    end
    diameter(str, :) = [featuresMSER.Axes(1), featuresMSER.Axes(2)];
    str = str + 1;
end
