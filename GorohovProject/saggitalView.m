%% Start with a clean slate
clear       %no variables
close all   %no figures
clc         %empty command window
imtool close all;

%% load MRI
load mri;
montage(D,map)
title('Horizontal Slices');

%% Mid saggital slice
M1 = D(:,64,:,:);
size(M1)
M2 = squeeze(M1);
size(M2)
%% show middle sagittal slice
figure
imshow(M2,map);
title('Sagittal - Raw data');