function moveitZ(h,A)
%MOVEITZ   Move a graphical object in 3-D in Z direction.
% Wang Gang, dnrwg@nus.edu.sg

% This code is modified from Anders Brun's moveit2 function.
% The original code can be found @ http://www.mathworks.com/matlabcentral/fileexchange/23304-moveit2-move-a-graphical-object-with-the-mouse

% Unpack gui object
gui = get(gcf,'UserData');

% Make a fresh figure window
set(h,'ButtonDownFcn',@startmovit);

% Store gui object
set(gcf,'UserData',{gui;A});


function startmovit(src,evnt)
% Unpack gui object
temp = get(gcf,'UserData');
gui=temp{1};
A=temp{2};

% Remove mouse pointer
set(gcf,'PointerShapeCData',nan(16,16));
set(gcf,'Pointer','custom');

% Set callbacks
gui.currenthandle = src;
thisfig = gcbf();
set(thisfig,'WindowButtonMotionFcn',@movit);
set(thisfig,'WindowButtonUpFcn',@stopmovit);

% Store starting point of the object
gui.startpoint = get(gca,'CurrentPoint');
set(gui.currenthandle,'UserData',{get(gui.currenthandle,'XData') get(gui.currenthandle,'YData') get(gui.currenthandle,'ZData')});

% Store gui object
set(gcf,'UserData',{gui;A});



function movit(src,evnt)
% Unpack gui object
temp = get(gcf,'UserData');
gui=temp{1};
A=temp{2};


try
if isequal(gui.startpoint,[])
    return
end
catch
end

% Do "smart" positioning of the object, relative to starting point...
pos = get(gca,'CurrentPoint')-gui.startpoint;
XYData = get(gui.currenthandle,'UserData');

set(gui.currenthandle,'ZData',XYData{3} + pos(1,3));
ZData=XYData{3};
ZData=max(ZData(:));
ZData=uint8(ZData + pos(1,3));
if ZData>size(A,3)
    ZData=size(A,3);
end
ImgCData=squeeze(A(:,:,ZData));

set(gui.currenthandle,'CData',ImgCData);
%set(gui.currenthandle,'CDataMapping','direct');


slicehandle=get(gui.currenthandle);
drawnow;

% Store gui object
set(gcf,'UserData',{gui;A});

function stopmovit(src,evnt)

% Clean up the evidence ...

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

