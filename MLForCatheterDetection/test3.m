openArea = 15;
ax = 'short'; % 'long1', 'long2'
img = GetImage(I, 48, ax);
[level, EM] = GetThresholdingLevel(img, 'ImprovedOtsu');
BW = imbinarize(img, level);
BW = bwareaopen(BW, openArea);
se = strel('ball', 3, 0, 0);
BW1 = imopen(BW, se);
% imshowpair(BW, BW1, 'montage');
% imshow(img, 'InitialMagnification', 'fit');
imshow(BW1, 'InitialMagnification', 'fit');