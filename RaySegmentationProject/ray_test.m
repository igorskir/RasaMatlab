clc;
clear;
%% Reading MRI images
functionPath = fullfile(pwd, 'Medical Image Reader and Viewer'); 
addpath(functionPath); 
fileFolder = fullfile(pwd, 'Gorohov'); 
CT = readImages(fileFolder); 
%---- end reading
%% Rays test
%preprocessing section: normalizaion and Gauss filter
image3D=permute(CT.volumes,[2 1 3]);
image3D=uint8(imnorm(image3D,'norm255'));
p = twomodegauss(0.1, 0.05, 0.3, 0.05, 1, 0.07, 0.002);
image3D=histeq(image3D,p);
 %setting  parameters...
points=[];
startPoint=[288, 165, 51];
dTetta=pi/8;
dFi=pi/8;
%processing
points=EmitRays(double(image3D),points,startPoint,dTetta,dFi,[]);
% marking up 3D Image
for i=1:size(points,1)
    image3D=putCross(image3D,[points(i,1) points(i,2) points(i,3)]);
end
implay(image3D,5);
imtool(squeeze(image3D(:,:,startPoint(3))));