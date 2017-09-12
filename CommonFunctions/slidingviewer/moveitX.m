function moveitX(h,A)
%MOVEITX   Move a graphical object in 3-D in X direction.
% Wang Gang, dnrwg@nus.edu.sg

% This code is modified from Anders Brun's moveit2 function.
% The original code can be found @ http://www.mathworks.com/matlabcentral/fileexchange/23304-moveit2-move-a-graphical-object-with-the-mouse


gui = get(gcf,'UserData');
set(h,'ButtonDownFcn',@startmovit);
set(gcf,'UserData',{gui;A});


function startmovit(src,evnt)
temp = get(gcf,'UserData');
gui=temp{1};
A=temp{2};

set(gcf,'PointerShapeCData',nan(16,16));
set(gcf,'Pointer','custom');

gui.currenthandle = src;
thisfig = gcbf();
set(thisfig,'WindowButtonMotionFcn',@movit);
set(thisfig,'WindowButtonUpFcn',@stopmovit);

gui.startpoint = get(gca,'CurrentPoint');
set(gui.currenthandle,'UserData',{get(gui.currenthandle,'XData') get(gui.currenthandle,'YData') get(gui.currenthandle,'ZData')});

set(gcf,'UserData',{gui;A});



function movit(src,evnt)
temp = get(gcf,'UserData');
gui=temp{1};
A=temp{2};


try
if isequal(gui.startpoint,[])
    return
end
catch
end

XYData = get(gui.currenthandle,'UserData');
 
pos = get(gca,'CurrentPoint')-gui.startpoint; 


set(gui.currenthandle,'XData',XYData{1} + pos(1,1));
XData=XYData{1};
XData=max(XData(:));
XData=round(XData + pos(1,1));
if XData > size(A,2)
    XData=size(A,2);
end
ImgCData=squeeze(A(:,XData,:)); 
set(gui.currenthandle,'CData',ImgCData);

slicehandle=get(gui.currenthandle);
drawnow;

set(gcf,'UserData',{gui;A});

function stopmovit(src,evnt)

thisfig = gcbf();

temp = get(gcf,'UserData');
gui=temp{1};
A=temp{2};

set(gcf,'Pointer','arrow');
set(thisfig,'WindowButtonUpFcn','');
set(thisfig,'WindowButtonMotionFcn','');
drawnow;
set(gui.currenthandle,'UserData','');
set(gcf,'UserData',{gui;A});

