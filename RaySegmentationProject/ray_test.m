%% Start with a clean slate
clear       %no variables
close all   %no figures
clc         %empty command window
imtool close all;
%% Reading MRI images
functionPath = fullfile(pwd, 'Medical Image Reader and Viewer'); 
addpath(functionPath); 
fileFolder = fullfile(pwd, 'Gorohov'); 
CT = readImages(fileFolder); 
%---- end reading

%% Data initialization
%preprocessing section: normalizaion and Gauss filter
image3D=permute(CT.volumes,[2 1 3]);
%image3D=uint8(imnorm(image3D,'norm255'));

%% Histogram modification
p = twomodegauss(0.1, 0.05, 0.3, 0.05, 1, 0.07, 0.002);
%image3D=histeq(image3D,p);

%% Setting  parameters
points=[];
startPoint=[245, 215, 51];
dTetta=pi/2;
dFi=pi/16;

%% Ray generation and processing
imtool(squeeze(image3D(:,:,startPoint(3))));
points=EmitRays(double(image3D),points,startPoint,dTetta,dFi);

%% Marking up 3D Image
for i=1:size(points,1)
    image3D=putCross(image3D,[points(i,1) points(i,2) points(i,3)]);
end

%% Show the results
%implay(image3D,2);
figure;
imtool(squeeze(image3D(:,:,startPoint(3))));
%K = convhulln(points(:,1:3));
%h=trisurf(K,points(:,1),points(:,2),points(:,3));
%set(h,'FaceColor','none');
