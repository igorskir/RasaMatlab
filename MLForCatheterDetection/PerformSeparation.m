% cd('Test watershed');
% load('BW16.mat');
% cd ..\
BW = bwareaopen(BW, 15, 8);
imshow(BW);
isVisual = 1;
scrSz = get(0, 'Screensize');
%%
D1 = -bwdist(~BW, 'euclidean');
D2 = -bwdist(~BW, 'cityblock');
D3 = -bwdist(~BW, 'chessboard');
D4 = -bwdist(~BW, 'quasi-Euclidean');

RGB1 = repmat(rescale(D1), [1 1 3]);
RGB2 = repmat(rescale(D2), [1 1 3]);
RGB3 = repmat(rescale(D3), [1 1 3]);
RGB4 = repmat(rescale(D4), [1 1 3]);

if isVisual == 1
    hFig = figure;
    subplot(2,2,1), imshow(RGB1), title('Euclidean')
    hold on, imcontour(D1)
    subplot(2,2,2), imshow(RGB2), title('City block')
    hold on, imcontour(D2)
    subplot(2,2,3), imshow(RGB3), title('Chessboard')
    hold on, imcontour(D3)
    subplot(2,2,4), imshow(RGB4), title('Quasi-Euclidean')
    hold on, imcontour(D4)
    set(gcf, 'Position', [scrSz(3), 0, scrSz(3), scrSz(4)],...
    'Color', 'w', 'numbertitle', 'off'); % FOR THE SECOND DISPLAY ONLY
end
%%
D1(~BW) = Inf;
D2(~BW) = Inf;
D3(~BW) = Inf;
D4(~BW) = Inf;
L1 = watershed(D1);
L1(~BW) = 0;
L2 = watershed(D2);
L2(~BW) = 0;
L3 = watershed(D3);
L3(~BW) = 0;
L4 = watershed(D4);
L4(~BW) = 0;

if isVisual == 1
    rgb1 = label2rgb(L1,'jet',[.5 .5 .5]);
    rgb2 = label2rgb(L2,'jet',[.5 .5 .5]);
    rgb3 = label2rgb(L3,'jet',[.5 .5 .5]);
    rgb4 = label2rgb(L4,'jet',[.5 .5 .5]);
    hFig = figure;
    subplot(2,2,1), imshowpair(rgb1, BW, 'montage'), title('Euclidean')
    subplot(2,2,2), imshowpair(rgb2, BW, 'montage'), title('City block')
    subplot(2,2,3), imshowpair(rgb3, BW, 'montage'), title('Chessboard')
    subplot(2,2,4), imshowpair(rgb4, BW, 'montage'), title('Quasi-Euclidean')
    set(gcf, 'Position', [scrSz(3), 0, scrSz(3), scrSz(4)],...
    'Color', 'w', 'numbertitle', 'off'); % FOR THE SECOND DISPLAY ONLY
end
%%
se = strel('disk', 2);
BW2 = logical(L1);
BWm = imopen(BW2, se);
CC = bwconncomp(BWm);
L = labelmatrix(CC);
cmap = colormap(brewermap([],'Accent'));
RGB = label2rgb(L, cmap, 'black', 'shuffle');
imshowpair(BW, RGB, 'montage');
set(gcf, 'Position', [scrSz(3), 0, scrSz(3), scrSz(4)],...
'Color', 'w', 'numbertitle', 'off'); % FOR THE SECOND DISPLAY ONLY
%%
if isVisual == 1
    stats = regionprops(L1,{'BoundingBox', 'Area', 'Extrema'});
    imshow(BW, 'InitialMagnification', 'fit');
    numWatershed = numel(stats);
    hold on
    for count = 1:numWatershed
        rectangle('Position', stats(count).BoundingBox, 'EdgeColor','c');
        posX = stats(count).Extrema(1,1) - 4;
        posY = stats(count).Extrema(1,2) - 3;
    end
    set(gcf, 'Position', scrSz, 'Color', 'w');
    hold off
end