% Otsu
level_otsu = graythresh(img);
BWadapt1 = imbinarize(img,level_otsu);
% Adaptive
T = adaptthresh(img,1,'ForegroundPolarity','dark');
BWadapt2 = imbinarize(img,T);
% Isodata
level_iso = isodata(img);
BWadapt3 = imbinarize(img, level_iso);
% Entropy th
BWadapt4 = im2bw_ent(img);

figure;
subplot(1,4,1);
imshow(BWadapt1);
addTitle('Otsu');

subplot(1,4,2);
imshow(BWadapt2, 'InitialMagnification', 'fit');
addTitle('Adaptive');

subplot(1,4,3);
imshow(BWadapt3);
addTitle('Isodata');

subplot(1,4,4);
imshow(BWadapt4);
addTitle('Entropy');