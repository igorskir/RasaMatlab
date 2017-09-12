I=zeros(176,176,208);
J1=I;
load 'J.mat';
J(:,:,20)=zeros(176,176,1);
J(:,:,69:72)=zeros(176,176,4);
J1(:,:,40:139)=J;
J=J1;
dataS = smooth3(J,'box',5);
dataS=J;
patch(isosurface(dataS,0.5),... 
'FaceColor',[0.745098 0.745098 0.745098],... 
'EdgeColor','none'); 
% isonormals(dataS,p1); 
%view(3); 
%axis vis3d tight 
view(156,30)
axis tight
axis([0 176 0 176 0 208])
daspect([1,1,1])
camlight right 
colormap('summer'); 
lighting gouraud

%PATCH_3Darray(dataS,1:h,1:w,1:n);
hold on;
grid on;
xlabel('X');
ylabel('Y');
zlabel('Z');
% PATCH_3Darray(permute(J,[2 1 3]));

I=zeros(176,176,208);

for i=1:208
    [m map]=imread(['input1/15/' num2str(i) '.bmp']);
    I(:,:,i)=m;
end;
I=I/255;
%I=implay(I);
level = graythresh(I);
for i=1:208
    I(:,:,i)=imbinarize(squeeze(I(:,:,i)),level);
end;
dataS = smooth3(I,'box',5);
patch(isosurface(dataS,0.9),... 
'FaceColor',[1 0.627451 0.478431],... 
'EdgeColor','none'); 