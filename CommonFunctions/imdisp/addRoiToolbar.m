function addRoiToolbar(fig)
%addRoiToolbar(fig)
%Add a toolbar for creating ROIs such as Line, Ellipse, Rectangle, Polygon and Freehand
%Right-click the ROIs to open a context menu, and you will see more
%functions there such as histogram, x-y plot, etc.
%Example:
%
%   load mri;
%   imshow(D(:,:,14))
%   addRoiToolbar;
%

%Yi Sui
%suiy02@gmail.com
%Apr 3rd, 2013

if nargin < 1
    fig=gcf;
end
figure(fig)

load addRoiToolbar.mat %icons
uipushtool('CData',icon{1},'TooltipString','Line','ClickedCallback',@roi_line_ClickedCallback);
uipushtool('CData',icon{2},'TooltipString','Ellipse','ClickedCallback',@roi_ellipse_ClickedCallback);
uipushtool('CData',icon{3},'TooltipString','Rectangle','ClickedCallback',@roi_rect_ClickedCallback);
uipushtool('CData',icon{4},'TooltipString','Polygon','ClickedCallback',@roi_poly_ClickedCallback);
uipushtool('CData',icon{5},'TooltipString','Freehand','ClickedCallback',@roi_freehand_ClickedCallback);
uipushtool('CData',icon{6},'TooltipString','Paste','ClickedCallback',@roi_paste_ClickedCallback);

% --------------------------------------------------------------------
function roi_line_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to roi_line (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

h=imdistline;
chld = get(h,'Children');
hcm=get(chld(end),'UIContextmenu');
item = uimenu(hcm, 'Label', 'X-Y Plot', 'Callback',@(varargin) roi_line_x_y_plot(h) );
uimenu(hcm, 'Label', 'Copy', 'Callback',@(varargin) roi_copy(h) );

% item = uimenu(hcm, 'Label', 'Delete', 'Callback',@(varargin) delete(h) );
addNewPositionCallback(h,@(varargin) roi_line_new_pos(h));
hline = findobj(h,'-regexp','Tag',' line');
btnfcn = get(hline,'ButtonDownFcn');

for k=1:numel(hline)
    
    set(hline(k), 'ButtonDownFcn',@(varargin) roi_line_btn_down_fcn(h,btnfcn{k}));
end

% --------------------------------------------------------------------
function roi_line_btn_down_fcn(h,fn)

mousebtn = get(gcf,'selectionType');
switch mousebtn
    case 'normal'
        roi_line_new_pos(h)
end
feval(fn)
% --------------------------------------------------------------------
function roi_line_new_pos(h)
ud=get(h,'UserData');

if isfield(ud,'xyplotwin') && ishandle(ud.xyplotwin) && strcmp(get(ud.xyplotwin,'type'),'figure')
    roi_line_x_y_plot(h)
end

% --------------------------------------------------------------------
function roi_line_x_y_plot (varargin)
% try
h=varargin{1};
img=getImage(h);
mask=createMask(h);
ud=get(h,'UserData');
if ~isfield (ud,'xyplotwin')
    ud.xyplotwin=figure;
    set(h,'Userdata',ud);
else
    figure(ud.xyplotwin);
end
plot(img(mask));


% --------------------------------------------------------------------
function roi_ellipse_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to roi_ellipse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h=imellipse(gca);
roi_init(h);

% --------------------------------------------------------------------
function roi_rect_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to roi_rect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

h=imrect;
roi_init(h);

% --------------------------------------------------------------------
function roi_poly_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to roi_poly (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)doc

h=impoly;
roi_init(h)


% --------------------------------------------------------------------
function roi_freehand_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to roi_freehand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h=imfreehand;
roi_init(h);

% --------------------------------------------------------------------

function roi_init(h)

chld = get(h,'Children');
hcm=get(chld,'UIContextmenu');
hcm= unique(cell2mat(hcm));
for k=1:length(hcm)
    uimenu(hcm(k), 'Label', 'Histogram', 'Callback',@(varargin) roi_hist(h) );
    uimenu(hcm(k), 'Label', 'Delete', 'Callback',@(varargin) roi_delete(h) );
    uimenu(hcm(k), 'Label', 'Copy', 'Callback',@(varargin) roi_copy(h) );
    uimenu(hcm(k), 'Label', 'Show/Hide Stats', 'Callback',@(varargin) show_hide_stats(h) );
end
setColor(h,'r')
pos=getPosition(h);
ud.htext = text(pos(1),pos(2),'','Tag','roi_text');
draggable(ud.htext);
ud.name = 'Unnamed';
% set(ud.htext,'ButtonDownFcn',@txtBDF)

set(h,'Userdata',ud);
ch = findobj(h,'type','line');
set(ch,'linewidth',1,'MarkerSize',3)

hpatch = findobj(h,'Tag','patch');
btnfcn = get(hpatch,'ButtonDownFcn');
set(hpatch, 'ButtonDownFcn',@(varargin) roi_btn_down_fcn(h,btnfcn));
roi_stat(h)
addNewPositionCallback(h,@(varargin) roi_new_pos(h));
%----------------------------------------------------
% function txtBDF(h)
% h

%-----------------------------------
function show_hide_stats(h) 
ud=get(h,'UserData');

if strcmpi( get(ud.htext,'Visible') , 'on')
    set(ud.htext,'Visible','off');
else
    set(ud.htext,'Visible','on');
end
    
%---------------------------------------------------------------------
function roi_copy(h)
assignin('base','roi_clipboard',h)
%----------------------------------
function roi_paste_ClickedCallback(hObject, eventdata, handles)
a=evalin('base','roi_clipboard');
pos = a.getPosition;
if strcmp(class(a),'imdistline')
    h = imdistline(gca,pos(:,1),pos(:,2));
    
elseif strcmp(class(a),'imfreehand')
    h = impoly(gca,pos);

else
    
    h = eval([class(a) '(gca, pos)']);
end
roi_init(h)

% --------------------------------------------------------------------
function roi_new_pos(h)
roi_stat(h)
ud=get(h,'UserData');
if isfield(ud,'histwin') && ishandle(ud.histwin) && strcmp(get(ud.histwin,'type'),'figure')
    roi_hist(h)
end
if isfield(ud,'xyplotwin') && ishandle(ud.xyplotwin) && strcmp(get(ud.xyplotwin,'type'),'figure')
    
    roi_line_x_y_plot(h)
    
end

% --------------------------------------------------------------------
function roi_btn_down_fcn(h,fn)
ud=get(h,'UserData');
obj=findobj('Tag','roi_text');
set(obj,'EdgeColor','none');
set(ud.htext,'EdgeColor',[1,0,0],'LineWidth',2);

mousebtn = get(gcf,'selectionType');
if isfield(ud,'histwin') && ishandle(ud.histwin)...
        && strcmp(get(ud.histwin,'type'),'figure')...
        && strcmp(mousebtn,'normal')
    roi_hist(h)
end
if strcmpi(mousebtn,'open')
    name=inputdlg('ROI name:','',1,{ud.name});
    if ~isempty(name)
    ud.name = name{1};
    set(h,'UserData',ud);
    roi_stat(h);
    end
end
feval(fn)

% --------------------------------------------------------------------
function roi_hist(varargin)
h=varargin{1};

ud=get(h,'UserData');
if ~isfield (ud,'histwin') || ~ishandle(ud.histwin)
    ud.histwin=figure;
    set(h,'Userdata',ud);
    
    h7 = uicontrol(...
        'Parent',ud.histwin,...
        'Units','normalized',...
        'BackgroundColor',[0.9 0.9 0.9],...
        'Callback',@(hObject,eventdata) histwin_callback(hObject,eventdata,h),...
        'Position',[0.18 0.01 0.65 0.04],...
        'String',{  'Slider' },...
        'Style','slider',...
        'CreateFcn',[],...
        'Tag','slider1');
    
else
    figure(ud.histwin)
end

hobj=findobj(ud.histwin,'tag','slider1');
histwin_callback(hobj,[],h)

% --------------------------------------------------------------------
function histwin_callback(hObject,eventdata,h)

N=5+get(hObject,'Value').*200;
% N
draw_hist(h,N)

function draw_hist(h,N)
if nargin<2
    N=100;
end

img=getImage(h);
mask=createMask(h);
hist(double(img(mask)),N);

% --------------------------------------------------------------------
function roi_stat(h)

img=getImage(h);
mask = createMask(h);
v=img(mask);
v=v(~isnan(v) & ~isinf(v));
v=double(v);
v_mean = mean(v);
v_std =std(v);
v_max = max(v);
v_min = min(v);

ud=get(h,'UserData');
ud.mean = v_mean;
ud.std = v_std;
ud.max = v_max;
ud.min = v_min;
ud.pix = sum(mask(:));
set(h,'UserData',ud);
set(ud.htext,...
    'String',sprintf('%s\nmean: %g, std: %g\nmin: %g, max:%g\nArea: %g pix (%g%%)',...
    ud.name,v_mean,v_std,v_min,v_max,sum(mask(:)),100*mean(mask(:))),...
    'BackgroundColor',[.8 .8 .8] );

% --------------------------------------------------------------------
function roi_delete(h)
ud=get(h,'UserData');
delete(ud.htext);
delete(h)

% --------------------------------------------------------------------
function img = getImage(h)
hax = get(h,'Parent');
himage=findobj(hax,'Type','image');
img = get(himage,'Cdata');