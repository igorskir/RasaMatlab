function slidingviewer(V)
%
% This function provides users, especially those doing medical image research, an
% interactive tool to browse through cross-sectional image scans (e.g. MRI and CT 
% images) from three orthogonal views (X,Y,Z) at their different levels.
% Users can click and hold on one image slice and move it to see the effect. 
% 
%%%% Demo %%%%
% load mri;
% D=double(squeeze(D));
% slidingviewer(D);
%%%%%%%%%%%%%% 
% 
%
%
%
% Author: Wang Gang, Department of Diagnostic Radiology,
% National University of Singapore
% Feedbacks are welcomed, please email to: dnrwg@nus.edu.sg
    figure;
    
    [x,y,z] = meshgrid(1:size(V,2),1:size(V,1),1:size(V,3));
     
	hslices=slice(x,y,z,V,size(V,2)/2,size(V,1)/2,size(V,3)/2);% slices by default @ middle pos.
    colormap(gray);
    xlabel('x');ylabel('y');zlabel('z');
    xlim([1 size(V,1)*4]); 
    ylim([1 size(V,2)*4]);
    zlim([1 size(V,3)*4]); 
    
    hXslice=hslices(1);
    hYslice=hslices(2);
    hZslice=hslices(3);
    
    set(hXslice,'EdgeColor','None','Tag','SliceX');
    set(hYslice,'EdgeColor','None','Tag','SliceY');
    set(hZslice,'EdgeColor','None','Tag','SliceZ');
    
    
    moveitX(hXslice,V);
    moveitY(hYslice,V);
    moveitZ(hZslice,V);
    
    
	axis tight;	
	axis vis3d;	
    set(gca,'ZDir','reverse');
    set(gca,'YDir','reverse');
      
end