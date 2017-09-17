clc;
clear;
%% Reading MRI images
functionPath = fullfile(pwd, 'Medical Image Reader and Viewer'); 
addpath(functionPath); 
fileFolder = fullfile(pwd, 'Gorohov'); 
CT = readImages(fileFolder); 

%---- end reading
%% Rays test

image3D=permute(CT.volumes,[2 1 3]);
image3D=uint8(imnorm(image3D,'norm255'));

EmitRays();
%p = twomodegauss(0.1, 0.05, 0.3, 0.05, 1, 0.07, 0.002); 
%image2D=histeq(image2D,p);


%image2D=squeeze(image3D(:,:,z0));
%image2D=squeeze(image3D(y0,:,:));% sagital
%imtool(image2D);
%imshow(squeeze(image3D(:,:,z0)));
%ray=Ray(image3D,x0,y0,z0,fi,tetta,delta);
%figure;
%plot(ray);
%figure;
%plot(ray);



% x=round(x0+dx*i);
% y=round(y0+dy*i);
% z=round(z0+dz*i);
% image2D=putCross(squeeze(image3D(:,:,z0)),x,y);
% imtool(image2D);