function show3D(data, isovalue)

% if islogical(data)
%     data = double(data);   
% end

sz = size(data);
surface = isosurface(data,isovalue);

scrSz = get(0, 'Screensize');
figure('Position', [scrSz(1), scrSz(2), scrSz(3)/2, scrSz(4)]);
addTitle('3D mesh');
hiso = patch('Vertices',surface.vertices,...
           'Faces',surface.faces,...
           'FaceColor',[1,.75,.65],...
           'FaceAlpha', 0.0,...
           'EdgeColor','red',...
           'EdgeAlpha', 0.15);
view(156,30)
axis tight
axis([0 sz(1) 0 sz(2) 0 sz(3)])
daspect([1,1,1])
set(gcf,'Renderer','opengl'); 
lighting none
set(gcf, 'Color', 'w');

figure('Position', [scrSz(3)/2, scrSz(2), scrSz(3)/2, scrSz(4)]);
addTitle('3D model');
hiso = patch('Vertices',surface.vertices,...
           'Faces',surface.faces,...
           'FaceColor',[1,.75,.65],...
           'FaceAlpha', 1.0,...
           'EdgeColor','black',...
           'EdgeAlpha', 0.0);
view(156,30) 
axis tight
axis([0 sz(1) 0 sz(2) 0 sz(3)])
daspect([1,1,1])
set(gcf,'Renderer','opengl'); 
lighting phong
material shiny
% lighting gouraud
camlight('headlight')
set(hiso,'SpecularColorReflectance',0,'SpecularExponent',50,'AmbientStrength',0.3)
set(gcf, 'Color', 'w');
end
