function show3Dcathinheart(data1, data2, isovalue1, isovalue2)

% if islogical(data)
%     data = double(data);   
% end

sz = size(data1);
surface1 = isosurface(data1,isovalue1);
surface2 = isosurface(data2,isovalue2);
scrSz = get(0, 'Screensize');

figure('Position', scrSz);
addTitle('3D model');
hiso1 = patch('Vertices',surface1.vertices,...
           'Faces',surface1.faces,...
           'FaceColor',[1,.75,.65],...
           'FaceAlpha', 1.0,...
           'EdgeColor','black',...
           'EdgeAlpha', 0.0);
% % view(156,30) 
% % axis tight
% % axis([0 sz(1) 0 sz(2) 0 sz(3)])
% % daspect([1,1,1])
% % set(gcf,'Renderer','opengl'); 
% % lighting phong
% % material shiny
% % % lighting gouraud
% % camlight('headlight')
% % set(hiso,'SpecularColorReflectance',0,'SpecularExponent',50,'AmbientStrength',0.3)
% % set(gcf, 'Color', 'w');

hiso2 = patch('Vertices',surface2.vertices,...
           'Faces',surface2.faces,...
           'FaceColor',[0.839216 0.839216 0.839216],...
           'FaceAlpha', 1.0,...
           'EdgeColor','black',...
           'EdgeAlpha', 0.0);
%0.839216 0.839216 0.839216
view(107,43) 
axis tight
axis([0 sz(1) 0 sz(2) 0 sz(3)])
daspect([1,1,1])
set(gcf,'Renderer','opengl');
lighting phong
material shiny
% lighting gouraud
camlight('headlight')
set(hiso1,'SpecularColorReflectance',0,'SpecularExponent',50,'AmbientStrength',0.3)
set(hiso2,'SpecularColorReflectance',0,'SpecularExponent',50,'AmbientStrength',0.3)
set(gcf, 'Color', 'w');

end
