% [imgCrop, rect] = imcrop(img);
k = 1.3;
level_not_crop = k*graythresh(img);
level_crop = graythresh(imgCrop);
BW_not_crop = imbinarize(imgCrop, level_not_crop);
BW_crop = imbinarize(imgCrop, level_crop);
imshowpair(BW_crop, BW_not_crop, 'montage');
%%
[T, V] = kittler(imgCrop);
BWk = imbinarize(imgCrop, T/255);
%%
figure
subplot(2,2,1);
imshow(imgCrop);
title('Image');

subplot(2,2,2);
imshow(BWk);
title('BWk');

subplot(2,2,3);
imshow(BW_crop);
title('BW_crop');

subplot(2,2,4);
imshow(BW_not_crop);
title('BW_not_crop');



