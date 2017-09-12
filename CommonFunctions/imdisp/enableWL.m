function enableWL(fh)
%Adjust Window/Level by mouse
if nargin<1
	fh=gcf;
end
% G=get(fh,'userdata');
G.oldWBMFcn = get(fh,'WindowButtonMotionFcn');
setappdata(fh,'enableWL_var',G);
set(fh,'WindowButtonDownFcn',@WBDFcn);
set(fh,'WindowButtonUpFcn',@WBUFcn);


function WBDFcn(varargin)
fh=varargin{1};
if ismember(get(fh,'SelectionType'),{'alt','extend'})
    G=getappdata(fh,'enableWL_var');
%     G.oldWBMFcn = get(fh,'WindowButtonMotionFcn');
    set(fh, 'WindowButtonMotionFcn',@AdjWL);
    G.initpnt=get(gca,'currentpoint');
    G.initClim = get(gca,'Clim');
    pt=get(gca,'CurrentPoint');
    try delete(G.txtbox);end
    G.txtbox=text(pt(1,1),pt(1,2),...
                  sprintf('W/L=[%.2g, %.2g]',G.initClim),...
                  'VerticalAlignment','bottom',...
                  'background','w');
%     uistack(gca,'top');
%     set(G.txtbox,'background','w');

    setappdata(fh,'enableWL_var',G);

    
end
    
function WBUFcn(varargin)
fh=varargin{1};
if ~strcmp(get(gcf,'SelectionType'),'normal')
G=getappdata(fh,'enableWL_var');
try delete(G.txtbox);end
set(fh,'WindowButtonMotionFcn',G.oldWBMFcn);
end


function AdjWL(varargin)
fh=varargin{1};
G=getappdata(fh,'enableWL_var');
G.cp=get(gca,'currentpoint');
G.x=G.cp(1,1);
G.y=G.cp(1,2);
G.xinit = G.initpnt(1,1);
G.yinit = G.initpnt(1,2);
G.dx = G.x-G.xinit;
G.dy = G.y-G.yinit;
G.clim = G.initClim+G.initClim(2).*[G.dx G.dy]./128;
pt=get(gca,'CurrentPoint');% 
set(G.txtbox,'position',[pt(1,1) pt(1,2)]);
set(G.txtbox,'string',sprintf('W/L=[%.2g, %.2g]',G.clim));
try
    switch get(fh,'SelectionType')
        case 'extend' % Mid-button, shft+left button,
        set(findobj(fh,'Type','axes'),'Clim',G.clim);
        case 'alt' %right-click,ctrl+left button,
        set(gca,'Clim',G.clim);
    end;

catch err
%     err.message
end;

