function saggitalviewer(V)

    
    [x,y,z] = meshgrid(1:size(V,2),1:size(V,1),1:size(V,3));
     
	hslices=slice(x,y,z,V,size(V,2)/2,size(V,1)/2,size(V,3)/2);% slices by default @ middle pos.
    colormap(gray);
    xlabel('x');ylabel('y');zlabel('z');
    xlim([1 size(V,1)*4]); 
    ylim([1 size(V,2)*4]);
    zlim([1 size(V,3)*4]); 
    
    %hXslice=hslices(1);
    hYslice=hslices(2);
    %hZslice=hslices(3);
    
    %set(hXslice,'EdgeColor','None','Tag','SliceX');
    set(hYslice,'EdgeColor','None','Tag','SliceY');
    %set(hZslice,'EdgeColor','None','Tag','SliceZ');
      
	axis tight;	
	axis vis3d;	
    set(gca,'ZDir','reverse');
    set(gca,'YDir','reverse');
      
end