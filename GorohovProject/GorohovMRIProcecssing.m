%% Start with a clean slate
clear       %no variables
close all   %no figures
clc         %empty command window
imtool close all;

addpath('D:\Projects\MatLabProjects\Common\slidingviewer');

%% Data Access
%----------------------------------------------------------------------

%filename convention used in image series (nobkpt)
% prefix = 'Gorohov\IMG0';
% fnum = 001:119;
% ext = '.dcm';
% 
% %first filename in series (nobkpt)
% fname = [prefix num2str(fnum(1)) ext];

fileFolder = fullfile(pwd, 'Gorohov');
files = dir(fullfile(fileFolder, '*dcm'));
fileNames = {files.name};

%examine file header (nobkpt)
info = dicominfo(fullfile(fileFolder, fileNames{1}));

%read one file to get size
I = dicomread(fullfile(fileFolder,fileNames{1}));
classI = class(I);
clear I;

%extract size info from metadata (nobkpt)
voxel_size = [info.PixelSpacing; info.SliceThickness];
numImages = length(fileNames);

%% Read slice images; populate 3-D matrix

D = zeros(info.Rows, info.Columns, numImages, classI);

%read slice images; populate XYZ matrix
hWaitBar = waitbar(0,'Reading DICOM files');
for i=numImages:-1:1
    fname = fullfile(fileFolder, fileNames{i});
    D(:,:,i) = uint16(dicomread(fname));
    waitbar((numImages-i)/numImages)
end
delete(hWaitBar)
whos D

%% Create Montage
figure
montage(reshape(uint16(D), [size(D,1), size(D,2), 1, numImages]),...
    'DisplayRange',[]);
%set(gca, 'clim', [0, 100]);
%drawnow;
shg; % Show most recent graph window


%% Exploring image data using Image Viewer GUI tool
im = squeeze(D(:,:,57));
max_level = double(max(im(:)));
imt = imtool(im, [0,max_level]);

%%% Exploring image using imdisp
%addpath('D:\Projects\MatLabProjects\Common\imdisp');
%%imdisp(D(:,:,20:40));

%% Custom display - image data
fig1 = figure;
fig1.Name = 'Coronal slice';
imshow(im,[0 max_level])
title('Coronal Slice #57')

colorbar     %add intensity legend 
colormap jet %change colormap


% %% Display Isosurface
% limits = [NaN NaN NaN NaN NaN 40];
% DS = squeeze(D);
% 
% [x, y, z, DS] = subvolume(DS, limits); 
% [fo,vo] = isosurface(x,y,z,DS,5);               % isosurface for the outside of the volume
% [fe,ve,ce] = isocaps(x,y,z,DS,5);               % isocaps for the end caps of the volume
% 
% %% Figurecreation
% figure
% p1 = patch('Faces', fo, 'Vertices', vo);       % draw the outside of the volume
% p1.FaceColor = 'red';
% p1.EdgeColor = 'none';
% 
% p2 = patch('Faces', fe, 'Vertices', ve, ...    % draw the end caps of the volume
%    'FaceVertexCData', ce);
% p2.FaceColor = 'interp';
% p2.EdgeColor = 'none';
% 
% view(-40,24)
% daspect([1 1 0.3])                             % set the axes aspect ratio
% colormap(gray(100))
% box on
% 
% camlight(40,40)                                % create two lights 
% camlight(-20,-10)
% lighting gouraud
%% Use Slidingviewer

slidingviewer(double(D));
%saggitalviewer(double(D));
%% Get Saggital plane
%T3 = tform = affine2d([Sx 0 0; 0 Sy 0; 0 0 Sz; 68.5 0 -14]);
Sx = -2.5;
Sy = 1;
Sz = 0.5;
tform = maketform('affine',[-2.5 0 0; 0 1 0; 0 0 0.5; 68.5 0 -14]);
R3 = makeresampler({'cubic','nearest','nearest'},'fill');
S = tformarray(D,tform,R3,[4 1 2],[1 2 4],[size(D,1) size(D,2) numImages],[],0);

%% Montage saggital

S2 = padarray(S(:,:,:,1),[2 0 0 0],0,'both');
% figure, montage(S2,[])
% title('Sagittal Slices');

montage(reshape(uint16(S2), [size(S2,1), size(S2,2), 1, numImages]),'DisplayRange',[]);
set(gca, 'clim', [0, 100]);
drawnow;
shg; % Show most recent graph window

%% Image binarisation
mriAdjust = D;
mriAdjust(280:end,:,:) = 0;
imshow(mriAdjust(:,:,57),[]);
bw = logical(mriAdjust);
imshow(bw(:,:,57),[]);

%% Image custing
A = imnorm(D,'norm255'); % нормализовали, на выходе double с макс знач 255
A = uint8(A); % изменили тип данных на uint8

%imtool(A(:,:,50),[]);
C = A(250,:,:);
imtool(A(250,:,:),[]);
%% 3D Visualisation section

%3D visualization (doc: contourslice, isosurface & isocap)
%docsearch('visualizing mri data')

% addpath('D:\Projects\MatLabProjects\Common\sliceomatic');
% sliceomatic(double(D));
% %ref: submission #780 @ www.mathworks.com/matlabcentral (nobkpt)
% hSlico1 = gcf;
% daspect(1./voxel_size)
% movegui('northwest')

%% intensity distribution also useful (more custom graphics)
%max_level = double(max(D(:))); 
my_map = jet(max_level);
fig2 = figure; 

%intensity distribution - top 2/3 (nobkpt)
subplot(3,1,1:2)
hist(double(im(:)),max_level)
axis([0 max_level 0 900])
title('Distribution')
