%% Start with a clean slate
clear       %no variables
close all   %no figures
clc         %empty command window
imtool close all;
% import square segmentation library
addpath('SqSegAlgo');
%% Reading MRI images
functionPath = fullfile(pwd, 'Medical Image Reader and Viewer'); 
addpath(functionPath); 
%functionPath = fullfile(pwd, 'Snake');  
%addpath(functionPath);
%functionPath = fullfile(pwd, 'BSpline'); 
%addpath(functionPath);
fileFolder = fullfile(pwd, 'Gorohov'); 
CT = readImages(fileFolder); 
%---- end reading

%% Data initialization
%preprocessing section: normalizaion and Gauss filter
image3D=permute(CT.volumes,[2 1 3]);
image3D=uint8(imnorm(image3D,'norm255'));

%% Histogram modification
p = twomodegauss(0.1, 0.05, 0.3, 0.05, 1, 0.07, 0.002);
image3D=histeq(image3D,p);

%% Setting  parameters
z=30; % slice number
image2D=squeeze(image3D(:,:,z));
%image2D=squeeze(sol_yxzt(:,:,5,12));
%addpath('testImages');
%image2D=rgb2gray(imread('brain2.jpg'));

imshow(image2D,[]);
coords=ginput(1);
startPoint= [coords(2) coords(1)];
% SqSeg options
options.threshold=78;
options.level=48;
options.edgeSize=2;
options.minEdgeSize=1; 

%SqSeg execution
tic();
splines=SquareSegmentation(image2D,startPoint,options);
toc();

figure;
imshow(image2D,[]);
externalContour=true;
for sp=1:numel(splines)
    splineContour=splines(sp); 
    t=0:0.7:splineContour.breaks(end);
    splinePoints=fnval(splineContour,t);
    x=[];
    y=[];
    for i=2:size(splinePoints,2)
        l=line([splinePoints(2,i-1) splinePoints(2,i)], [splinePoints(1,i-1) splinePoints(1,i)]);
        x(end+1)=splinePoints(1,i);
        y(end+1)=splinePoints(2,i);
        set(l, {'color','LineWidth'}, {'red',2});
    end
    l=line([splinePoints(2,i) splinePoints(2,1)], [splinePoints(1,i) splinePoints(1,1)]);
    set(l, {'color','LineWidth'}, {'red',2});
    x(end+1)=splinePoints(1,1);
    y(end+1)=splinePoints(2,1);
    if(externalContour)
        mask=poly2mask(x,y,size(image2D,2),size(image2D,1))';
        externalContour=false;
    else
        mask=mask-poly2mask(x,y,size(image2D,2),size(image2D,1))';
    end
end
 figure;
 imshow(mask);

