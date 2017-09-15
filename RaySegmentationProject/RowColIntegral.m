clc;
clear;

%load 'MRI data/sol_yxzt_pat18.mat';
%A=uint8(squeeze(sol_yxzt(:,:,5,12)));

path = 'MRI data\Gorohov\03310959\'; 
image3D=dicomread([path 'data (' num2str(1) ')']);
[height, width]=size(image3D);
slices=95;
image3D=zeros(height,width,slices);
for i=1:slices
    image3D(:,:,i)=dicomread([path 'data (' num2str(i) ')']);  
end
image3D=uint8(imnorm(image3D,'norm255'));
margin=100;
image3D=image3D(margin:height-margin,margin:width-margin,:);

%implay(image3D,10);
%[M,N] = size(A); 
%ColSum = zeros(M,1); 
%RowSum = zeros(N,1); 
%A = double(A); 
%fileFolder='Gorohov\03310959';
%%
fileFolder = fullfile(pwd, 'Gorohov'); 
CT = readImages(fileFolder); 
% PET volume is already aligned to scanner coordinate axes 


VolumeViewer3D(CT); 
%%
% tic(); 
% for i = 1:M 
% for j = 1:N 
%     ColSum(j) = ColSum(j) + A(i,j); 
%     RowSum(i) = RowSum(i) + A(i,j); 
% end 
% end 
% toc(); 

%Region Growing Testing

%J=regiongrowing(double(A)/255,M/2,N/2,0.1);
%imshow(A);
%figure;
%imshow(J);



%3D Region Growing Testing

%A=uint8(squeeze(sol_yxzt(:,:,:,12)));



%H=[1 1 0; 1 1 0; 0 0 0];
%A=zeros(3,3,3);
%A(:,:,1)=H;
%A(:,:,2)=H;
%A(:,:,3)=H;
%A
% [M N K]=size(A);
% 
% tic();
%     %J=regiongrowing3D(double(A)/255,140,125,4,0.7);
%     [P J]=regionGrowing3D_slow(double(A)/255,[140 125 4]);
% toc()

%implay(J,3);
%implay(A,3);
%for i=1:K
%    montage(squeeze(J(:,:,i)),squeeze(A(:,:,i)));
%end;
% I=image3D;
% level = graythresh(I);
% for i=slices
%     I(:,:,i)=imbinarize(squeeze(I(:,:,i)),level);
% end;
% dataS=I;
% patch(isosurface(dataS,0.4),... 
% 'FaceColor',[0.745098 0.745098 0.745098],... 
% 'EdgeColor','none'); 
% % isonormals(dataS,p1); 
% %view(3); 
% %axis vis3d tight 
% view(156,30)
% axis tight
% %axis([0 176 0 176 0 208])
% daspect([1,1,1])
% camlight right 
% colormap('summer'); 
% lighting gouraud
% 
% %PATCH_3Darray(dataS,1:h,1:w,1:n);
% hold on;
% grid on;
% xlabel('X');
% ylabel('Y');
% zlabel('Z');

