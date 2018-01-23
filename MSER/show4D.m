function show4D(input, framerate)
pauseTime = 1/framerate; 
scrSz = get(0, 'Screensize');
sz = size(input);
% figure('Position', [scrSz(3)/2, scrSz(2), scrSz(3)/2, scrSz(4)]);
figure;
set(gcf, 'Position', scrSz, 'Color', 'w');
addTitle('3D model');
surface = getSurfaces(input, 0.075);
envFlag = 0;
stopFlag = 1;
pause(0.5);
while stopFlag == 1
    for i = 1:sz(4)
    hiso = patch('Vertices',surface(i).vertices,...
               'Faces',surface(i).faces,...
               'FaceColor',[1,.75,.65],...
               'FaceAlpha', 1.0,...
               'EdgeColor','black',...
               'EdgeAlpha', 0.0);
        if envFlag ~= 1      
            camlight;
            view(115,50) % for the initial reconstruction
%             view(156,30) % for the final reconstruction
            axis tight
            axis([0 sz(1) 0 sz(2) 0 sz(3)])
            daspect([1,1,1])
            camproj orthographic
            set(gcf,'Renderer','opengl');
            % lighting phong
            material shiny
            hiso.FaceLighting = 'gouraud';
            hiso.BackFaceLighting = 'reverselit';
            hiso.AmbientStrength = 0.3;
            hiso.DiffuseStrength = 0.65;
            hiso.SpecularStrength = 1;
            hiso.SpecularExponent = 50;
            hiso.SpecularColorReflectance = 0;
            set(gcf, 'Color', 'w');

        end   
    envFlag = 1;
    str = sprintf('3D model\n Timeframe: %d', i);
    addTitle(str);
    pause(pauseTime);
    delete(hiso)
    end
end
end

