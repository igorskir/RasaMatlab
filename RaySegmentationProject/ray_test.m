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


tetta=pi/2; %initial values
fi=-pi/2;
delta=1;
x0=255;
y0=255;
z0=51;
imshow(squeeze(image3D(:,:,z0)));
ray=Ray(image3D,x0,y0,z0,fi,tetta,delta);
%figure;
%plot(ray);
ray=smooth(ray,5);
%figure;
%plot(ray);

stride=3;
ray2=circshift(ray,stride);
rayDiff=abs(ray-ray2);
rayDiff(1:stride)=0;
rayDiff(numel(rayDiff)-stride:end)=0;
i=1;
while(rayDiff(i)<25 && i<numel(rayDiff)) 
    i=i+1;
end
dx=delta*sin(tetta)*cos(fi); % projection step on each dimension 
dy=-delta*sin(tetta)*sin(fi);%
dz=delta*cos(tetta);%
x=round(x0+dx*i);
y=round(y0+dy*i);
z=round(z0+dz*i);
image2D=putCross(squeeze(image3D(:,:,z0)),x,y);
imtool(image2D);