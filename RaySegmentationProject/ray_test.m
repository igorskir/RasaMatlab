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
%fi=-pi/2;
delta=1;
x0=225;
y0=225;
z0=51;
width=2;
level=1.15;
%image2D=squeeze(image3D(:,:,z0));
image2D=squeeze(image3D(y0,:,:));% sagital

p = twomodegauss(0.1, 0.05, 0.3, 0.05, 1, 0.07, 0.002); 
image2D=histeq(image2D,p);
fi=0;
for tetta=-pi/2:pi/16:pi/2
    ray=Ray(image3D, x0,y0,z0,fi,tetta,delta);
    borderIndex=GetBorder(ray,delta,width,level);
    if(borderIndex==-1)
           continue;
    end
    dx=delta*sin(tetta)*cos(fi); % projection step on each dimension 
    dy=-delta*sin(tetta)*sin(fi);%
    dz=delta*cos(tetta);%
    x=round(x0+dx*borderIndex);
    y=round(y0+dy*borderIndex);
    z=round(z0+dz*borderIndex);
    image2D=putCross(image2D,x,y);
end
image2D=putCross(image2D,x0,y0);
imtool(image2D);
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