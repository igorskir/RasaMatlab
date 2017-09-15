clc;
clear;
%% Reading MRI images
fileFolder = fullfile(pwd, 'Gorohov'); 
CT = readImages(fileFolder); 

%---- end reading
%% RG test
image3D=permute(CT.volumes,[2 1 3]);
%implay(uint8(image3D));
segmentedMap=regiongrowing3D(double(image3D)/255,205,242,51,0.2);
implay(segmentedMap);
%%  Visualization
patch(isosurface(double(segmentedMap),0.2),... 
'FaceColor',[0.745098 0.745098 0.745098],... 
'EdgeColor','none'); 
view(156,30)
axis tight
%axis([0 176 0 176 0 208])
daspect([1,1,1])
camlight right 
colormap('summer'); 
lighting flat