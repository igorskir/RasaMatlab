
function varargout = VolumeViewer3D(varargin)

%-------------------------------------------------------------------------%
%                                                                         %
% code by JDS, 2015                                                       %
%                                                                         %
% Possible inputs are                                                     %
%   {DIRECTORY}                                                           %
%   {VOLUME}                                                              %
%   {VOLUME, voxDimX, voxDimY, voxDimZ}                                   %
%   {STRUCTURE} (FROM 'readImages' FUNCTION)                              %
%   {STRUCTURE, STRUCTURE} for fused volume viewing                       %
%   {..., [specific application tagword]} (OPTIONAL)                      %
%       -'align'  aligns volume to scanner coordinates                    %
%       -'trim'  trim volume to shared space (only for fused viewing)     %
%       -'full'  use full size image matrices (only for fused viewing)    %
%       -'getCoord'  returns user selected voxel location                 %
%       -'drawROI' (or 'drawVOI')  returns a user drawn mask              %
%                                                                         %
%-------------------------------------------------------------------------%

addpath(genpath(cd))

warning('off','MATLAB:warn_r14_stucture_assignment')

global handles
handles = struct;

%fused image display is only supported for two structure inputs
handles.display.images.fusedImages = false;
%voxel location will be displayed if using structure(s) or reading from file
handles.display.text.showLocs = false;
%calculating mask?
handles.display.images.drawROI = false;
handles.mask.draw = false; %toggles if drawing

makeInactive()
drawnow

%PARSE INPUTS
inStruct = parseInputs(varargin);
if isempty(inStruct) || (isnumeric(inStruct) && inStruct==-1)
    if nargout>0, varargout = cell(1,nargout); end
    return
elseif isstruct(inStruct)
    if (~isfield(inStruct,'vol') && ~isfield(inStruct,'vol_1')) || ...
            (isfield(inStruct,'vol') && isempty(inStruct.vol))
        disp('Incomplete or unrecognized inputs')
        if nargout>0, varargout = cell(1,nargout); end
        return
    end
end

if isfield(inStruct,'vol_1') && isfield(inStruct,'vol_2')
    handles.display.images.fusedImages = true;
    if isfield(inStruct,'units_1')
        handles.display.images.units_1 = inStruct.units_1;
    else
        handles.display.images.units_1 = '';
    end
    if isfield(inStruct,'units_2')
        handles.display.images.units_2 = inStruct.units_2;
    else
        handles.display.images.units_2 = '';
    end
    handles.volume.modality_1 = inStruct.modality_1;
    handles.volume.modality_2 = inStruct.modality_2;
else
    if isfield(inStruct,'units')
        handles.display.images.units = inStruct.units;
    else
        handles.display.images.units = '';
    end
    handles.volume.modality = inStruct.modality;
end
if isfield(inStruct,'bedRange1') && isfield(inStruct,'bedRange2') && isfield(inStruct,'bedRange3')
    handles.volume.bedRange1 = inStruct.bedRange1;
    handles.volume.bedRange2 = inStruct.bedRange2;
    handles.volume.bedRange3 = inStruct.bedRange3;
elseif isfield(inStruct,'voxDim1') && isfield(inStruct,'voxDim2') && isfield(inStruct,'voxDim3')
    [handles.volume.bedRange1,handles.volume.bedRange2,handles.volume.bedRange3] = ...
        getVolumeBedRanges(inStruct);
end
if isfield(inStruct,'patientOrientation')
    handles.volume.patientOrientation = inStruct.patientOrientation;
end
if isfield(inStruct,'drawROI')
    handles.display.images.drawROI = true;
end

%the following increases the voxel value range for display purposes
%sometimes necessary for non-quantitative images
if ~handles.display.images.fusedImages
    volMax = max(inStruct.vol(:));
    volMin = min(inStruct.vol(:));
    if volMin==0 && volMax==0
        close(findobj('Tag','3DviewerFig'))
        disp('Volume has no non-zero elements')
        if nargout>0, varargout = cell(1,nargout); end
        return
    end
    handles.display.images.sclFctr = 1;
    if max(inStruct.vol(:))<=0
        relScl = -1;
    else
        relScl = 1;
    end
    inStruct.vol = inStruct.vol * relScl;
    volMax = volMax * relScl;
    volMin = volMin * relScl;    
    while max([volMin volMax]) < 100
        inStruct.vol = inStruct.vol * 100;
        handles.display.images.sclFctr = handles.display.images.sclFctr*100;
        volMax = volMax * 100;
        volMin = volMin * 100;
    end    
    inStruct.vol = inStruct.vol * relScl;
    handles.display.images.volMax = volMax * relScl;
    handles.display.images.volMin = volMin * relScl;
else
    volMax_1 = max(inStruct.vol_1(:));
    volMin_1 = min(inStruct.vol_1(:));
    if volMin_1==0 && volMax_1==0
        close(findobj('Tag','3DviewerFig'))
        disp('Input 1 has no non-zero elements')
        if nargout>0, varargout = cell(1,nargout); end
        return
    end
    volMax_2 = max(inStruct.vol_2(:));
    volMin_2 = min(inStruct.vol_2(:));
    if volMin_2==0 && volMax_2==0
        close(findobj('Tag','3DviewerFig'))
        disp('Input 2 has no non-zero elements')
        if nargout>0, varargout = cell(1,nargout); end
        return
    end
    
    handles.display.images.sclFctr_1 = 1;
    if max(inStruct.vol_1(:))<=0
        relScl_1 = -1;
    else
        relScl_1 = 1;
    end
    inStruct.vol_1 = inStruct.vol_1 * relScl_1;
    volMax_1 = volMax_1 * relScl_1;
    volMin_1 = volMin_1 * relScl_1;    
    while max([volMin_1 volMax_1]) < 100
        inStruct.vol_1 = inStruct.vol_1 * 100;
        handles.display.images.sclFctr_1 = handles.display.images.sclFctr_1*100;
        volMax_1 = volMax_1 * 100;
        volMin_1 = volMin_1 * 100;
    end    
    inStruct.vol_1 = inStruct.vol_1 * relScl_1;
    handles.display.images.volMax_1 = volMax_1 * relScl_1;
    handles.display.images.volMin_1 = volMin_1 * relScl_1;
    
    handles.display.images.sclFctr_2 = 1;
    if max(inStruct.vol_2(:))<=0
        relScl_2 = -1;
    else
        relScl_2 = 1;
    end
    inStruct.vol_2 = inStruct.vol_2 * relScl_1;
    volMax_2 = volMax_2 * relScl_2;
    volMin_2 = volMin_2 * relScl_2;    
    while max([volMin_2 volMax_2]) < 100
        inStruct.vol_2 = inStruct.vol_2 * 100;
        handles.display.images.sclFctr_2 = handles.display.images.sclFctr_2*100;
        volMax_2 = volMax_2 * 100;
        volMin_2 = volMin_2 * 100;
    end    
    inStruct.vol_2 = inStruct.vol_2 * relScl_1;
    handles.display.images.volMax_2 = volMax_2 * relScl_2;
    handles.display.images.volMin_2 = volMin_2 * relScl_2;
end

if isfield(inStruct,'filePath')
    handles.volume.filePath = inStruct.filePath;
else
    handles.volume.filePath = '{unknown}';
end
%INITIALIZE FIGURE, AXES, SLIDERS, AND TEXT BOXES
initialSetup(inStruct);
figSz = get(handles.figure.f1,'Position');
figSz = figSz(3:4);
%initialize colormaps after figure (only needed for fused images)
if handles.display.images.fusedImages
    test1 = get(findobj('Tag','CMPanel_1'),'Children');
    test2 = get(findobj('Tag','CMPanel_2'),'Children');
    switch lower(inStruct.cm_1)
        case 'gray'
            if strcmpi(handles.volume.modality_1,'pet')
                handles.display.images.cm_1 = 1-gray;
            else
                handles.display.images.cm_1 = gray;
            end
            set(test1(strcmpi(get(test1,'String'),'gray')),'Value',1)
        case 'hot'
            handles.display.images.cm_1 = hotMetal;
            set(test1(strcmpi(get(test1,'String'),'hot')),'Value',1)
        case 'pet'
            handles.display.images.cm_1 = PET;
            set(test1(strcmpi(get(test1,'String'),'PET')),'Value',1)
    end
    switch lower(inStruct.cm_2)
        case 'gray'
            if strcmpi(handles.volume.modality_2,'pet')
                handles.display.images.cm_2 = 1-gray;
            else
                handles.display.images.cm_2 = gray;
            end
            set(test2(strcmpi(get(test2,'String'),'gray')),'Value',1)
        case 'hot'
            handles.display.images.cm_2 = hotMetal;
            set(test2(strcmpi(get(test2,'String'),'hot')),'Value',1)
        case 'pet'
            handles.display.images.cm_2 = PET;
            set(test2(strcmpi(get(test2,'String'),'PET')),'Value',1)
    end
end

if isfield(inStruct,'IOP')
    handles.volume.IOP = inStruct.IOP;
    if all(round(inStruct.IOP*1e4)/1e4==[1;0;0;0;1;0])
        title(handles.axes.a1,{'Transaxial';''},'FontSize',0.01*figSz(1));
        title(handles.axes.a2,{'Coronal';''},'FontSize',0.01*figSz(1));
        title(handles.axes.a3,{'Sagittal';''},'FontSize',0.01*figSz(1));
    else
        title(handles.axes.a1,'','FontSize',0.01*figSz(1));
        title(handles.axes.a2,'','FontSize',0.01*figSz(1));
        title(handles.axes.a3,'','FontSize',0.01*figSz(1));
    end
else
    title(handles.axes.a1,'','FontSize',0.01*figSz(1));
    title(handles.axes.a2,'','FontSize',0.01*figSz(1));
    title(handles.axes.a3,'','FontSize',0.01*figSz(1));
end

makeInactive()

handles.display.images.fltr = 1;
handles.display.images.fltr = handles.display.images.fltr/sum(handles.display.images.fltr(:));

if isfield(inStruct,'drawROI')
    %VOLUME PROPERTIES
    loadVolume(inStruct);
    %INITIALIZE MASK, PATCH DISPLAY, AND SLICE INDECES
    handles.mask.mask = false(handles.volume.imSz1,handles.volume.imSz2,handles.volume.imSz3);
    axes(handles.axes.a1)
    handles.mask.a1patch = patch([0;0;0;0;0],[0;0;0;0;0],...
        'g','EdgeColor',[0.5 1 0],'FaceAlpha',0.1,'EdgeAlpha',0.8,'HitTest','off','Tag','maskPatch1');
    axes(handles.axes.a2)
    handles.mask.a2patch = patch([0;0;0;0;0],[0;0;0;0;0],...
        'g','EdgeColor',[0.5 1 0],'FaceAlpha',0.1,'EdgeAlpha',0.8,'HitTest','off','Tag','maskPatch2');
    axes(handles.axes.a3)
    handles.mask.a3patch = patch([0;0;0;0;0],[0;0;0;0;0],...
        'g','EdgeColor',[0.5 1 0],'FaceAlpha',0.1,'EdgeAlpha',0.8,'HitTest','off','Tag','maskPatch3');
    try
        set(handles.mask.a1patch,'PickableParts','none')
        set(handles.mask.a2patch,'PickableParts','none')
        set(handles.mask.a3patch,'PickableParts','none')
    catch
    end
    [handles.mask.a1Yidx,handles.mask.a1Xidx] = meshgrid(1:handles.volume.imSz2,1:handles.volume.imSz1);
    [handles.mask.a2Yidx,handles.mask.a2Xidx] = meshgrid(1:handles.volume.imSz3,1:handles.volume.imSz1);
    [handles.mask.a3Yidx,handles.mask.a3Xidx] = meshgrid(1:handles.volume.imSz3,1:handles.volume.imSz2);
    %DRAW BUTTON
    defineDrawButton()
    %DONE BUTTON
    defineDoneButton()
    %MAKE FIGURE RESIZEABLE
    makeFigScalable()
    %DEFINE ALL CALLBACKS
    makeActive()
    uiwait();
    
    varargout = cell(1,nargout);
    
    if isempty(findobj('Tag','3DviewerFig'))
        disp('There is no output')
        return;
    end
    close(handles.figure.f1)
    drawnow
    
    mask = handles.mask.mask;
    
    %transform mask into original space if needed
    if isfield(inStruct,'locX')
        %get original image volume for reference
        for n=1:numel(varargin)
            if isstruct(varargin{n}) && (isfield(varargin{n},'volume') || ...
                    isfield(varargin{n},'volumes') || isfield(varargin{n},'vol'))
                origVol = varargin{n};
                break
            end
        end
        if ~exist('origVol','var')
            %file name or directory input
            origVol = inStruct;
        end
        if isfield(origVol,'volume')
            siz = size(origVol.volume);
        elseif isfield(origVol,'volumes')
            siz = size(origVol.volumes);
        elseif isfield(origVol,'vol')
            siz = size(origVol.vol);
        end
        if ~(numel(size(inStruct.locX))==numel(siz) && all(size(inStruct.locX)==siz)) || ...
                ~all(inStruct.bedRange1==origVol.bedRange1) || ...
                ~all(inStruct.bedRange2==origVol.bedRange2) || ...
                ~all(inStruct.bedRange3==origVol.bedRange3)
            %transform voxel space
            inStruct.volumes = mask;
            temp = transformVoxelSpace(inStruct,origVol);
            mask = temp.volumes;
        end
    end
    varargout{1} = mask;
elseif isfield(inStruct,'getCoord')
    %VOLUME PROPERTIES
    loadVolume(inStruct);
    %DONE BUTTON
    defineDoneButton()
    %MAKE FIGURE RESIZEABLE
    makeFigScalable()
    %DEFINE ALL CALLBACKS
    makeActive()
    uiwait();
    
    varargout = cell(1,nargout);
        
    if isempty(findobj('Tag','3DviewerFig'))
        disp('There is no output')
        return;
    end
    close(handles.figure.f1)
    drawnow
    
    %transform coordinates into original space if needed
    if isfield(inStruct,'locX')
        %get original volume for reference
        for n=1:numel(varargin)
            if isstruct(varargin{n})
                origVol = varargin{n};
                break
            end
        end
        if ~exist('origVol','var')
            %file name or directory input
            origVol = inStruct;
        end
        if isfield(origVol,'volume')
            siz = size(origVol.volume);
        elseif isfield(origVol,'volumes')
            siz = size(origVol.volumes);
        elseif isfield(origVol,'vol')
            siz = size(origVol.vol);
        end
        if ~(numel(size(inStruct.locX))==numel(siz) && all(size(inStruct.locX)==siz)) || ...
                ~all(inStruct.bedRange1==origVol.bedRange1) || ...
                ~all(inStruct.bedRange2==origVol.bedRange2) || ...
                ~all(inStruct.bedRange3==origVol.bedRange3)
            %coordinates in real dimensions (mm)
            x = inStruct.locX(handles.display.images.x,handles.display.images.y,handles.display.images.z);
            y = inStruct.locY(handles.display.images.x,handles.display.images.y,handles.display.images.z);
            z = inStruct.locZ(handles.display.images.x,handles.display.images.y,handles.display.images.z);
            
            %check for coordinates outside of original coverage, accounting for original IOP
            xr = origVol.IOP(1);
            yr = origVol.IOP(2);
            zr = origVol.IOP(3);
            xc = origVol.IOP(4);
            yc = origVol.IOP(5);
            zc = origVol.IOP(6);
            xs = yr*zc-yc*zr;
            ys = zr*xc-zc*xr;
            zs = xr*yc-xc*yr;
            rot3D = [xr xc xs; yr yc ys; zr zc zs];
            temp = rot3D'*[origVol.locX(1,1,1); origVol.locY(1,1,1); origVol.locZ(1,1,1)];
            rotRange1(1) = temp(1);
            rotRange2(1) = temp(2);
            rotRange3(1) = temp(3);
            temp = rot3D'*[origVol.locX(end,1,1); origVol.locY(end,1,1); origVol.locZ(end,1,1)];
            rotRange1(2) = temp(1);
            temp = rot3D'*[origVol.locX(1,end,1); origVol.locY(1,end,1); origVol.locZ(1,end,1)];
            rotRange2(2) = temp(2);
            temp = rot3D'*[origVol.locX(1,1,end); origVol.locY(1,1,end); origVol.locZ(1,1,end)];
            rotRange3(2) = temp(3);
            rotPt = rot3D'*[x; y; z];
            voxDim = [abs(diff(origVol.bedRange1))/size(origVol.locX,1) ...
                abs(diff(origVol.bedRange2))/size(origVol.locX,2) ...
                abs(diff(origVol.bedRange3))/size(origVol.locX,3)];
            if rotPt(1)>max(rotRange1)+voxDim(1)/2 || rotPt(1)<min(rotRange1)-voxDim(1)/2 || ...
                    rotPt(2)>max(rotRange2)+voxDim(2)/2 || rotPt(2)<min(rotRange2)-voxDim(2)/2 || ...
                    rotPt(3)>max(rotRange3)+voxDim(3)/2 || rotPt(3)<min(rotRange3)-voxDim(3)/2
                disp('Selected location appears to be outside of original volume')
                varargout{1} = [];
                return
            end
            %find voxel in original volume closest to the new location
            if ~isfield(origVol,'locX')
                [origVol.locX,origVol.locY,origVol.locZ] = getVoxelLocs(origVol);
            end
            distVolSq = (origVol.locX-x).^2+(origVol.locY-y).^2+(origVol.locZ-z).^2;
            idx = find(distVolSq==min(distVolSq(:))); idx = idx(1);
            %convert index into coordinates in original volume
            siz = size(origVol.locX);
            k = [1 cumprod(siz(1:end-1))];
            for n = numel(siz):-1:1,
                vi = rem(idx-1, k(n)) + 1;
                vj = (idx - vi)/k(n) + 1;
                switch n
                    case 1
                        handles.display.images.x = double(vj);
                    case 2
                        handles.display.images.y = double(vj);
                    case 3
                        handles.display.images.z = double(vj);
                end
                idx = vi;
            end
            if numel(siz)<3, handles.display.images.z = 1; end
            if numel(siz)<2, handles.display.images.y = 1; end
        end
    end
    varargout{1} = [handles.display.images.x, handles.display.images.y, handles.display.images.z];
elseif isfield(inStruct,'vol') || handles.display.images.fusedImages
    %VOLUME
    loadVolume(inStruct);
    %MAKE FIGURE RESIZEABLE
    makeFigScalable()
    %DEFINE ALL CALLBACKS
    makeActive()
    %mitigate error by passing empty output
    if nargout>0
        disp('There is no output')
        varargout = cell(1,nargout);
    end
else
    error('Incorrect inputs')
end

end

function out = parseInputs(in)

allowedTags = {'align','trim','full','getcoord','getcoords','drawroi','drawvoi'};

out = struct; %default return is changed to -1 if no image(s) loaded

align = false;
fusedImages = false;

SUVunits = false;

%% check inputs
structIdx = []; %count number of structure inputs to see if they came from 'readImages' function
numlogIdx = []; %count number of numerical or logical structures (e.g. voxDim, mask input, etc.)
strIdx = []; %count number of string inputs (e.g. image directory, tag words, etc.)
rmIdx = []; %inputs to be removed

numIdx = 0; %keep track of extra inputs to accompany image data
for n=1:numel(in)
    rmFlag = false;
    
    if isstruct(in{n})
        %structure inputs
        structIdx = [structIdx n];
        if ~isfield(in{n},'volume') && ~isfield(in{n},'volumes')
            disp(['Image data is not defined in (structure) input ' num2str(n)])
            rmFlag = true;
        elseif ~isfield(in{n},'IOP')
            disp(['Image Orientation Patient is not defined in (structure) input ' num2str(n)])
            rmFlag = true;
        elseif ~isfield(in{n},'IPP')
            disp(['Image Position Patient is not defined in (structure) input ' num2str(n)])
            rmFlag = true;
        elseif ~isfield(in{n},'bedRange1') && ~isfield(in{n},'voxDim1')
            disp(['Cannot get voxel dimension 1 from (structure) input ' num2str(n)])
            rmFlag = true;
        elseif ~isfield(in{n},'bedRange2') && ~isfield(in{n},'voxDim2')
            disp(['Cannot get voxel dimension 2 from (structure) input ' num2str(n)])
            rmFlag = true;
        elseif ~isfield(in{n},'bedRange3') && ~isfield(in{n},'voxDim3')
            disp(['Cannot get voxel dimension 3 from (structure) input ' num2str(n)])
            rmFlag = true;
        end
        if rmFlag
            disp(['   ...so ignoring (structure) input ' num2str(n)])
        end
    elseif isnumeric(in{n}) || islogical(in{n})
        %numerical or logical inputs
        if isempty(numlogIdx) || numlogIdx(end)+numIdx<=n
            numlogIdx = [numlogIdx n];
            if numel(size(in{n}))==4 || numel(size(in{n}))==3 || ...
                    (numel(size(in{n}))==2 && numel(in{n})>1)
                %image input
                if (numel(in)>=n+1 && isnumeric(in{n+1}) && numel(in{n+1})==1) && ...
                        (numel(in)>=n+2 && isnumeric(in{n+2}) && numel(in{n+2})==1) && ...
                        (numel(in)>=n+3 && isnumeric(in{n+3}) && numel(in{n+3})==1)
                    %3D voxel dimensions found
                    numlogIdx = [numlogIdx n+1];
                    numlogIdx = [numlogIdx n+2];
                    numlogIdx = [numlogIdx n+3];
                    numIdx = numIdx+3;
                end
            else
                if islogical(in{n+numIdx})
                    disp(['The (logical) input ' num2str(n+numIdx) ' is not image or voxel data matching'])
                    disp(['   ...so ignoring (logical) input ' num2str(n)])
                else
                    disp(['The (numerical) input ' num2str(n+numIdx) ' is not image or voxel data'])
                    disp(['   ...so ignoring (numerical) input ' num2str(n)])
                end
                rmFlag = true;
            end
        else
            numIdx = numIdx-1;
        end
    elseif ischar(in{n})
        %string inputs
        origStr = in{n};
        in{n} = strtrim(in{n});
        strIdx = [strIdx n];
        switch lower(in{n})
            case (allowedTags)
                rmFlag = false;
            otherwise
                if ~(exist(in{n},'file')==2) && ~(exist(in{n},'file')==7)
                    disp(['''' origStr ''' is not a recognized file, directory, or flag'])
                    disp(['   ...so ignoring (string) input ' num2str(n)])
                    rmFlag = true;
                end
        end
    else
        disp(['Input type ''' class(in{n}) ''' is not allowed'])
        disp(['   ...so ignoring (' class(in{n}) ') input ' num2str(n)])
        rmFlag = true;
    end
    if rmFlag
        %index removed input
        rmIdx = [rmIdx n];
    end
end

%cleanup ignored inputs
in(rmIdx) = [];
structIdx(ismember(structIdx,rmIdx)) = [];
numlogIdx(ismember(numlogIdx,rmIdx)) = [];
strIdx(ismember(strIdx,rmIdx)) = [];
%then adjust indeces
for n=1:numel(rmIdx)
    structIdx(structIdx>=rmIdx(n)) = structIdx(structIdx>=rmIdx(n))-1;
    numlogIdx(numlogIdx>=rmIdx(n)) = numlogIdx(numlogIdx>=rmIdx(n))-1;
    strIdx(strIdx>=rmIdx(n)) = strIdx(strIdx>=rmIdx(n))-1;
    rmIdx(rmIdx>=rmIdx(n)) = rmIdx(rmIdx>=rmIdx(n))-1;
end
rmIdxOrig = rmIdx;

%% record and remove options, and perform checks
varargin = {};
chkRpt = {};
rmIdx = [];
if any(strcmpi(in,'align'))
    align = true;
    rmIdx = [rmIdx find(strcmpi(in,'align'))];
end
if any(strcmpi(in,'drawROI')) || any(strcmpi(in,'drawVOI'))
    if any(strcmpi(in,'getCoord')) || any(strcmpi(in,'getCoords'))
        chkRpt{numel(chkRpt)+1} = 'Ignoring ''getCoord'' flag in favor of drawing ROI mask';
        rmIdx = [rmIdx find(strcmpi(in,'getCoord') | strcmpi(in,'getCoords'))];
    end
    out.drawROI = true;
    rmIdx = [rmIdx find(strcmpi(in,'drawROI') | strcmpi(in,'drawVOI'))];
else
    if any(strcmpi(in,'getCoord')) || any(strcmpi(in,'getCoords'))
        out.getCoord = true;
        rmIdx = [rmIdx find(strcmpi(in,'getCoord') | strcmpi(in,'getCoords'))];
    end
end
if any(strcmpi(in,'trim'))
    varargin = [varargin, {'trim'}];
    rmIdx = [rmIdx find(strcmpi(in,'trim'))];
end
if any(strcmpi(in,'full'))
    varargin = [varargin, {'full'}];
    rmIdx = [rmIdx find(strcmpi(in,'full'))];
end
rmIdx = sort(unique(rmIdx));

%cleanup string inputs
in(rmIdx) = [];
structIdx(ismember(structIdx,rmIdx)) = [];
numlogIdx(ismember(numlogIdx,rmIdx)) = [];
strIdx(ismember(strIdx,rmIdx)) = [];
%then adjust indeces
for n=1:numel(rmIdx)
    structIdx(structIdx>=rmIdx(n)) = structIdx(structIdx>=rmIdx(n))-1;
    numlogIdx(numlogIdx>=rmIdx(n)) = numlogIdx(numlogIdx>=rmIdx(n))-1;
    strIdx(strIdx>=rmIdx(n)) = strIdx(strIdx>=rmIdx(n))-1;
    rmIdx(rmIdx>=rmIdx(n)) = rmIdx(rmIdx>=rmIdx(n))-1;
end

%first check if directory passed so can return a structure
if ~isempty(strIdx)
    disp('Reading data from all valid input files or directories')
    rmIdx = [];
    for n=1:numel(strIdx)
        try
            temp = readImages(in{strIdx(n)});
            if isempty(temp) || (isnumeric(temp) && temp==-1)
                disp(['There was a problem with input ''' in{strIdx(n)} ''''])        
                numlogIdx = sort([numlogIdx strIdx(n)]);
            else
                %check if multiple series returned
                if isfield(temp,'nSeries') && temp.nSeries>1
                    disp('Multiple series found in directory, using the first one')
                    fn = fieldnames(temp);
                    fn(strcmpi(fn,'nSeries')) = [];
                    temp = temp.(fn{1});
                end                
                structIdx = sort([structIdx strIdx(n)]);
            end
            %now replace this string input with the structure
            in{strIdx(n)} = temp;
            rmIdx = [rmIdx n];
        catch
            %character string input is not file or directory (but snuck by previous check)
            disp(['''' in{strIdx(n)} ''' is not a recognized file or directory'])
            disp('   ...so ignoring it')
        end
    end
    strIdx(rmIdx) = [];
    %cleanup remaining string inputs
    rmIdx = strIdx;
    in(rmIdx) = [];
    structIdx(ismember(structIdx,rmIdx)) = [];
    numlogIdx(ismember(numlogIdx,rmIdx)) = [];
    strIdx(ismember(strIdx,rmIdx)) = [];
    %then adjust indeces
    for n=1:numel(rmIdx)
        structIdx(structIdx>=rmIdx(n)) = structIdx(structIdx>=rmIdx(n))-1;
        numlogIdx(numlogIdx>=rmIdx(n)) = numlogIdx(numlogIdx>=rmIdx(n))-1;
        strIdx(strIdx>=rmIdx(n)) = strIdx(strIdx>=rmIdx(n))-1;
        rmIdx(rmIdx>=rmIdx(n)) = rmIdx(rmIdx>=rmIdx(n))-1;
    end
end

if numel(structIdx)>1
    fusedImages = true;
    SUVunits_1 = false;
    SUVunits_2 = false;
end
if numel(structIdx)>=1 && numel(numlogIdx)>0
    %this allows a structure with a matching mask volume to be input
    siz = [];
    if isfield(in{structIdx(1)},'volume')
        siz = size(in{structIdx(1)}.volume);
    elseif isfield(in{structIdx(1)},'volumes')
        siz = size(in{structIdx(1)}.volumes);
    end
    for n=1:numel(numlogIdx)
        temp = size(in{numlogIdx(n)});
        if numel(siz)==numel(temp) && all(siz==temp)
            %format into basic structure and replace numerical input
            volStruct = in{structIdx(1)};
            if isfield(volStruct,'units'), volStruct = rmfield(volStruct,'units'); end
            if isfield(volStruct,'modality'), volStruct = rmfield(volStruct,'modality'); end
            volStruct.volumes = in{numlogIdx(n)};
            structIdx = [structIdx numlogIdx(n)];
            in{numlogIdx(n)} = volStruct;
            numlogIdx(n) = [];
            fusedImages = true;
            SUVunits_1 = false;
            SUVunits_2 = false;
            break
        end
    end
end
if numel(structIdx)>2
    disp('Cannot fuse more than 2 volumes')
%     out = -1;
%     return
end

if align
    if ~isempty(structIdx)
        if isfield(in{structIdx(1)},'volume')
            siz = size(in{structIdx(1)}.volume);
        elseif isfield(in{structIdx(1)},'volumes')
            siz = size(in{structIdx(1)}.volumes);
        end
        if numel(siz)==2
            chkRpt{numel(chkRpt)+1} = 'Ignoring ''align'' flag for 2D slice input';
        else
            varargin = [varargin, {'align'}];
            align = true;
        end
    end
end
if isfield(out,'drawROI') && out.drawROI
    if fusedImages
        %check if both 2D and 3D inputs
        if isfield(in{structIdx(1)},'volume')
            siz1 = size(in{structIdx(1)}.volume);
        elseif isfield(in{structIdx(1)},'volumes')
            siz1 = size(in{structIdx(1)}.volumes);
        end
        if isfield(in{structIdx(2)},'volume')
            siz2 = size(in{structIdx(2)}.volume);
        elseif isfield(in{structIdx(2)},'volumes')
            siz2 = size(in{structIdx(2)}.volumes);
        end
        if numel(siz1)~=numel(siz2)
            %if 3D volume was first input, it should not be transformed to 2D space
            if numel(siz1)>numel(siz2)
                chkRpt{numel(chkRpt)+1} = 'Ignoring 2D slice input since defining 3D volume mask';
                rmIdx = [rmIdx structIdx(2)];
                fusedImages = false;
            end
        end
    end
end

%display unused input check report
for n=1:numel(chkRpt)
    disp(chkRpt{n})
end
%this condition does not cover every situation
if (numel(numlogIdx)>0 && numel(structIdx)>0) || ...
        numel(numlogIdx)==2 || numel(numlogIdx)==3 || numel(numlogIdx)>4
    if isempty(rmIdx)
        disp('Ignoring unused inputs')
    else
        disp('Ignoring other unused inputs')
    end
end

%arrange inputs to give priority to structures
in = in([sort(structIdx) sort(numlogIdx) sort(strIdx)]);
structIdx = 1:numel(structIdx);
numlogIdx = (1:numel(numlogIdx))+numel(structIdx);
strIdx = (1:numel(strIdx))+numel(structIdx)+numel(numlogIdx);

%no input
if isempty(in)
    in = readImages();
    if isempty(in) || (isnumeric(in) && in==-1)
        disp('There was a problem finding or reading the image files')
        out = -1;
        return
    end
    %check if multiple series returned
    if isfield(in,'nSeries') && in.nSeries>1
        disp('Multiple series found, using the first one')
        fn = fieldnames(in);
        fn(strcmpi(fn,'nSeries')) = [];
        in = in.(fn{1});
    end
    if numel(size(in.volumes))==4
        disp('Displaying last frame of dynamic series')
        in.volumes = in.volumes(:,:,:,end);
    end
    %convert structure to cell to match 'in'
    temp = cell(1);
    temp{1} = in;
    in = temp;
end

for n=1:numel(in)
    if isnumeric(in{n}) || islogical(in{n})
        siz = size(in{n});
        if siz(1)~=1 && siz(2)~=1 %must be slice or volume
            out.modality = 'VOL';
            switch numel(size(in{n}))
                case {2, 3} %slice or volume
                    out.vol = in{n};
                case 4
                    disp('Displaying last frame of dynamic series')
                    out.vol = in{n}(:,:,:,end);
            end
            if align
                disp('Ignoring ''align'' flag because no structure input')
            end
            if numel(in)>=n+3 && ...
                    isnumeric(in{n+1}) && isnumeric(in{n+2}) && isnumeric(in{n+3}) && ...
                    numel(in{n+1})==1 && numel(in{n+2})==1 && numel(in{n+3})==1
                out.FOV1 = single(in{n+1} * size(out.vol, 1));
                out.FOV2 = single(in{n+2} * size(out.vol, 2));
                out.FOV3 = single(in{n+3} * size(out.vol, 3));
            else
                %generic FOVs
                out.FOV1 = single(size(out.vol, 1));
                out.FOV2 = single(size(out.vol, 2));
                out.FOV3 = single(size(out.vol, 3));
            end
            return
        end
    elseif isstruct(in{n}) %ASSUME STRUCTURE FROM FUNCTION 'readImages'
        %check if multiple series in input structure
        if isfield(in{n},'nSeries') && in{n}.nSeries>1
            disp('Multiple series in input structure, using the first series')
            fn = fieldnames(in{n});
            fn(strcmpi(fn,'nSeries')) = [];
            in{n} = in{n}.(fn{1});
        end
        %convert to SUV if all necessary information exists
        if isfield(in{n},'injectedDose'), injDose = in{n}.injectedDose; end
        if isfield(in{n},'injectionTime'), injTime = in{n}.injectionTime; end
        if isfield(in{n},'injectionDate'), injDate = in{n}.injectionDate; end
        if isfield(in{n},'studyTime'), acqTime = in{n}.studyTime; end
        if isfield(in{n},'studyDate'), acqDate = in{n}.studyDate; end
        if isfield(in{n},'units'), actUnits = in{n}.units; end
        if isfield(in{n},'tracer'), tracer = in{n}.tracer; end
        if isfield(in{n},'patientWeight'), patWeight = in{n}.patientWeight; end
        try
            SUVconv = ...
                calculateSUV(injDose, injTime, injDate, acqTime, acqDate, patWeight, actUnits, tracer);
            if ~isnan(SUVconv) && ~isinf(SUVconv)
                SUVunits = true;
            end
        catch
        end
        if ~fusedImages
            out.modality = 'VOL';
            temp = parseInputStruct(in{n});
            
            if isfield(temp,'modality')
                out.modality = temp.modality;
                if strcmpi(out.modality,'PT'), out.modality = 'PET'; end
            else
                if isfield(temp,'gantryModel')
                    if strcmpi(temp.gantryModel,'1104') || strcmpi(temp.gantryModel(end-2:end),'mCT')
                        out.modality = 'CT';
                    elseif strcmpi(temp.gantryModel,'2008') || strcmpi(temp.gantryModel(end-2:end),'mMR')
                        out.modality = 'MR';
                    end
                end
                if isfield(temp,'units') && strcmpi(temp.units,'HU')
                    out.modality = 'CT';
                end
                if isfield(temp,'units') && strcmpi(temp.units,'1/CM')
                    out.modality = 'umap';
                end
                if isfield(temp,'tracer') || ...
                        (isfield(temp,'units') && strcmpi(temp.units,'BQML'))
                    out.modality = 'PET';
                end
            end
            %check for binary mask data
            if isfield(temp,'volume')
                if islogical(temp.volume)
                    out.modality = 'mask';
                end
            elseif isfield(temp,'volumes')
                if islogical(temp.volumes)
                    out.modality = 'mask';
                end
            end
            if align
                if ~all(round(temp.IOP*1e6)/1e6==[1;0;0;0;1;0]) || ...
                        temp.bedRange1(1)>temp.bedRange1(2) || ...
                        temp.bedRange2(1)>temp.bedRange2(2) || ...
                        temp.bedRange3(1)>temp.bedRange3(2)
                    disp(['   ' out.modality ' will be aligned to scanner coordinate axes'])
                    if isfield(out,'drawROI') && out.drawROI
                        disp(['   ...but output mask will be in voxel space of ' out.modality ' input, i.e. not ''aligned'' space'])
                    elseif isfield(out,'getCoord') && out.getCoord
                        disp(['   ...but output coordinates will be in voxel space of ' out.modality ' input, i.e. not ''aligned'' space'])
                    end
                    temp.volumes = temp.vol;
                    temp1 = alignVolume(temp);
                    disp('      aligning volume...')
                    if isfield(temp1,'volume')
                        temp.vol = temp1.volume;
                    elseif isfield(temp1,'volumes')
                        temp.vol = temp1.volumes;
                    end
                    temp.IOP = temp1.IOP;
                    temp.IPP = temp1.IPP;
                    temp.bedRange1 = single(temp1.bedRange1);
                    temp.bedRange2 = single(temp1.bedRange2);
                    temp.bedRange3 = single(temp1.bedRange3);
                    temp.FOV1 = single(abs(diff(temp1.bedRange1)));
                    temp.FOV2 = single(abs(diff(temp1.bedRange2)));
                    temp.FOV3 = single(abs(diff(temp1.bedRange3)));
                    temp.locX = single(temp1.locX);
                    temp.locY = single(temp1.locY);
                    temp.locZ = single(temp1.locZ);
                else
                    [bedRange1, bedRange2, bedRange3, locX, locY, locZ] = ...
                        getVolumeBedRanges(temp);
                    temp.bedRange1 = single(bedRange1);
                    temp.bedRange2 = single(bedRange2);
                    temp.bedRange3 = single(bedRange3);
                    temp.FOV1 = single(abs(diff(bedRange1)));
                    temp.FOV2 = single(abs(diff(bedRange2)));
                    temp.FOV3 = single(abs(diff(bedRange3)));
                    temp.locX = single(locX);
                    temp.locY = single(locY);
                    temp.locZ = single(locZ);
                end
            else
                [bedRange1, bedRange2, bedRange3, locX, locY, locZ] = ...
                    getVolumeBedRanges(temp);
                temp.bedRange1 = single(bedRange1);
                temp.bedRange2 = single(bedRange2);
                temp.bedRange3 = single(bedRange3);
                temp.FOV1 = single(abs(diff(bedRange1)));
                temp.FOV2 = single(abs(diff(bedRange2)));
                temp.FOV3 = single(abs(diff(bedRange3)));
                temp.locX = single(locX);
                temp.locY = single(locY);
                temp.locZ = single(locZ);
            end
            fn = fieldnames(temp);
            for m=1:numel(fn)
                out.(fn{m}) = temp.(fn{m});
            end
            if strcmpi(out.modality,'PT'), out.modality = 'PET'; end
            if SUVunits
                out.vol = out.vol * SUVconv;
                out.units = 'SUV';
            end
            if islogical(out.vol)% || all(out.vol(:)==0 | out.vol(:)==1)
                out.modality = 'mask';
                if isfield(out,'units')
                    out = rmfield(out,'units');
                end
            end
        else
            if exist('SUVconv','var')
                SUVconv_1 = SUVconv;
                SUVunits_1 = true;
            end
            struct1 = in{n};
            for m=n+1:numel(in)
                if isstruct(in{m})
                    %check if multiple series in input structure
                    if isfield(in{m},'nSeries') && in{m}.nSeries>1
                        disp('Multiple series in second input structure, using the first series')
                        fn = fieldnames(in{m});
                        fn(strcmpi(fn,'nSeries')) = [];
                        in{m} = in{m}.(fn{1});
                    end
                    struct2 = in{m};
                    if isfield(struct2,'injectedDose'), injDose_2 = struct2.injectedDose; end
                    if isfield(struct2,'injectionTime'), injTime_2 = struct2.injectionTime; end
                    if isfield(struct2,'injectionDate'), injDate_2 = struct2.injectionDate; end
                    if isfield(struct2,'studyTime'), acqTime_2 = struct2.studyTime; end
                    if isfield(struct2,'studyDate'), acqDate_2 = struct2.studyDate; end
                    if isfield(struct2,'units'), actUnits_2 = struct2.units; end
                    if isfield(struct2,'tracer'), tracer_2 = struct2.tracer; end
                    if isfield(struct2,'patientWeight'), patWeight_2 = struct2.patientWeight; end
                    try
                        SUVconv_2 = ...
                            calculateSUV(injDose_2, injTime_2, injDate_2, acqTime_2, acqDate_2, patWeight_2, actUnits_2, tracer_2);
                        if ~isnan(SUVconv_2) && ~isinf(SUVconv_2)
                            SUVunits_2 = true;
                        end
                    catch
                    end
                    break
                end
            end
            %use the combination of data to calculate SUV conversions (only patient weight)
            if isfield(struct1,'tracer')
                if ~SUVunits_1
                    try
                        SUVconv_1 = ...
                            calculateSUV(injDose, injTime, injDate, acqTime, acqDate, patWeight_2, actUnits, tracer);
                        if ~isnan(SUVconv_1) && ~isinf(SUVconv_1)
                            SUVunits_1 = true;
                        end
                    catch
                    end
                end
            end
            if isfield(struct2,'tracer')
                if ~SUVunits_2
                    try
                        SUVconv_2 = ...
                            calculateSUV(injDose_2, injTime_2, injDate_2, acqTime_2, acqDate_2, patWeight, actUnits_2, tracer_2);
                        if ~isnan(SUVconv_2) && ~isinf(SUVconv_2)
                            SUVunits_2 = true;
                        end
                    catch
                    end
                end
            end
            %image modalities
            out.modality_1 = 'VOL';
            out.modality_2 = 'VOL';
            out.cm_1 = 'gray';
            out.cm_2 = 'gray';
            if isfield(struct1,'modality')
                out.modality_1 = struct1.modality;
                if strcmpi(out.modality_1,'PT'), out.modality_1 = 'PET'; end
            else
                if isfield(struct1,'gantryModel')
                    if strcmpi(struct1.gantryModel,'1104') || strcmpi(struct1.gantryModel(end-2:end),'mCT')
                        out.modality_1 = 'CT';
                        out.cm_1 = 'gray';
                    elseif strcmpi(struct1.gantryModel,'2008') || strcmpi(struct1.gantryModel(end-2:end),'mMR')
                        out.modality_1 = 'MR';
                        out.cm_1 = 'gray';
                    end
                elseif isfield(struct1,'units') && strcmpi(struct1.units,'HU')
                    out.modality_1 = 'CT';
                    out.cm_1 = 'gray';
                end
            end
            if isfield(struct2,'modality')
                out.modality_2 = struct2.modality;
                if strcmpi(out.modality_2,'PT'), out.modality_2 = 'PET'; end
            else
                if isfield(struct2,'gantryModel')
                    if strcmpi(struct2.gantryModel,'1104') || strcmpi(struct2.gantryModel(end-2:end),'mCT')
                        out.modality_2 = 'CT';
                        out.cm_2 = 'gray';
                    elseif strcmpi(struct2.gantryModel,'2008') || strcmpi(struct2.gantryModel(end-2:end),'mMR')
                        out.modality_2 = 'MR';
                        out.cm_2 = 'gray';
                    end
                elseif isfield(struct2,'units') && strcmpi(struct2.units,'HU')
                    out.modality_2 = 'CT';
                    out.cm_2 = 'gray';
                end
            end
            if isfield(struct1,'tracer') || (isfield(out,'modality_1') && strcmpi(out.modality_1,'PET'))
                out.modality_1 = 'PET';
                out.cm_1 = 'PET';
            end
            if isfield(struct2,'tracer') || (isfield(out,'modality_2') && strcmpi(out.modality_2,'PET'))
                out.modality_2 = 'PET';
                out.cm_2 = 'PET';
            end
            if isfield(struct1,'units') && strcmpi(struct1.units,'1/CM')
                out.modality_1 = 'umap';
                out.cm_1 = 'gray';
            end
            if isfield(struct2,'units') && strcmpi(struct2.units,'1/CM')
                out.modality_2 = 'umap';
                out.cm_2 = 'gray';
            end
            %check for binary mask data
            if isfield(struct1,'volume')
                if islogical(struct1.volume)
                    out.modality_1 = 'mask';
                    out.cm_1 = 'gray';
                end
            elseif isfield(struct1,'volumes')
                if islogical(struct1.volumes)
                    out.modality_1 = 'mask';
                    out.cm_1 = 'gray';
                end
            end
            if isfield(struct2,'volume')
                if islogical(struct2.volume)
                    out.modality_2 = 'mask';
                    out.cm_2 = 'gray';
                end
            elseif isfield(struct2,'volumes')
                if islogical(struct2.volumes)
                    out.modality_2 = 'mask';
                    out.cm_2 = 'gray';
                end
            end
            %check for dynamic series
            if isfield(struct1,'volume')
                if size(struct1.volume,4)>1
                    disp(['Displaying last frame of dynamic ' out.modality_1 ' series'])
                    struct1.volume = struct1.volume(:,:,:,end);
                end
            elseif isfield(struct1,'volumes')
                if size(struct1.volumes,4)>1
                    disp(['Displaying last frame of dynamic ' out.modality_1 ' series'])
                    struct1.volumes = struct1.volumes(:,:,:,end);
                end
            end
            if isfield(struct2,'volume')
                if size(struct2.volume,4)>1
                    disp(['Displaying last frame of dynamic ' out.modality_2 ' series'])
                    struct2.volume = struct2.volume(:,:,:,end);
                end
            elseif isfield(struct2,'volumes')
                if size(struct2.volumes,4)>1
                    disp(['Displaying last frame of dynamic ' out.modality_2 ' series'])
                    struct2.volumes = struct2.volumes(:,:,:,end);
                end
            end
            if isfield(struct1,'IOP') && isfield(struct2,'IOP')
                if ~isfield(struct1,'locX') || ~isfield(struct1,'locY') || ~isfield(struct1,'locZ')
                    [struct1.locX,struct1.locY,struct1.locZ] = getVoxelLocs(struct1);
                end
                if ~isfield(struct2,'locX') || ~isfield(struct2,'locY') || ~isfield(struct2,'locZ')
                    [struct2.locX,struct2.locY,struct2.locZ] = getVoxelLocs(struct2);
                end                
            else
                disp('Cannot coregister inputs because image orientation axes are not defined')
                out = -1;
                return
            end
            
            %check if dimensions of inputs to determine coregistration method
            if isfield(in{structIdx(1)},'volume')
                siz1 = size(in{structIdx(1)}.volume);
            elseif isfield(in{structIdx(1)},'volumes')
                siz1 = size(in{structIdx(1)}.volumes);
            end
            if isfield(in{structIdx(2)},'volume')
                siz2 = size(in{structIdx(2)}.volume);
            elseif isfield(in{structIdx(2)},'volumes')
                siz2 = size(in{structIdx(2)}.volumes);
            end
            if ~isfield(out,'drawROI') && numel(siz1)==numel(siz2) %mutually coregister
                if align
                    if ~all(round(in{structIdx(1)}.IOP*1e6)/1e6==[1;0;0;0;1;0]) || ...
                        in{structIdx(1)}.bedRange1(1)>in{structIdx(1)}.bedRange1(2) || ...
                        in{structIdx(1)}.bedRange2(1)>in{structIdx(1)}.bedRange2(2) || ...
                        in{structIdx(1)}.bedRange3(1)>in{structIdx(1)}.bedRange3(2)
                        disp('   Volumes will be aligned to scanner coordinate axes')
                        if isfield(out,'getCoord') && out.getCoord
                            if strcmpi(out.modality_1,out.modality_2)
                                disp(['   ...but output coordinates will be in voxel space of first ' out.modality_1 ' input, i.e. not ''aligned'' space'])
                            else
                                disp(['   ...but output coordinates will be in voxel space of ' out.modality_1 ' input, i.e. not ''aligned'' space'])
                            end
                        end
                    end
                else
                    %check relative alignments
                    if ~all(round(struct2.IOP*1e6)/1e6==round(struct1.IOP*1e6)/1e6) || ...
                            ((struct2.bedRange1(1)>struct2.bedRange1(2) && struct1.bedRange1(1)<struct1.bedRange1(2)) || ...
                            (struct2.bedRange1(1)<struct2.bedRange1(2) && struct1.bedRange1(1)>struct1.bedRange1(2))) || ...
                            ((struct2.bedRange2(1)>struct2.bedRange2(2) && struct1.bedRange2(1)<struct1.bedRange2(2)) || ...
                            (struct2.bedRange2(1)<struct2.bedRange2(2) && struct1.bedRange2(1)>struct1.bedRange2(2))) || ...
                            ((struct2.bedRange3(1)>struct2.bedRange3(2) && struct1.bedRange3(1)<struct1.bedRange3(2)) || ...
                            (struct2.bedRange3(1)<struct2.bedRange3(2) && struct1.bedRange3(1)>struct1.bedRange3(2)))
                        %do not need to align here, it will be handled by coregistration functions
                        if strcmpi(out.modality_1,out.modality_2)
                            disp(['   Second ' out.modality_2 ' input will be aligned to orientation of first'])
                            if isfield(out,'getCoord') && out.getCoord
                                disp(['   ...but output coordinates will be in voxel space of first ' out.modality_1 ' input, i.e. not coregistered space'])
                            end
                        else
                            disp(['   ' out.modality_2 ' will be aligned to orientation of ' out.modality_1])
                            if isfield(out,'getCoord') && out.getCoord
                                disp(['   ...but output coordinates will be in voxel space of ' out.modality_1 ' input, i.e. not coregistered space'])
                            end
                        end
                    end
                    %both volumes will be aligned to IOP from first input
                end
                disp('      coregistering volumes...')
                temp = coregisterVolumes(struct1,struct2,varargin{:});
                if isempty(temp) || (isnumeric(temp) && temp==-1)
                    out = -1;
                    return
                end
                out.IOP = temp.IOP;
            else %transform one to the other
                if isfield(out,'drawROI')
                    %match voxel grid if drawing VOI mask
                    if ~align
                        if strcmpi(out.modality_1,out.modality_2)
                            disp(['   The second ' out.modality_2 ' input will be matched to the voxel space of the first for drawing mask VOI'])
                        else
                            disp(['   The ' out.modality_2 ' input will be matched to the voxel space of the ' out.modality_1 ' for drawing mask VOI'])
                        end
                    else
                        if isfield(in{structIdx(1)},'volume')
                            siz = size(in{structIdx(1)}.volume);
                        elseif isfield(in{structIdx(1)},'volumes')
                            siz = size(in{structIdx(1)}.volumes);
                        end
                        if numel(siz)>2
                            if strcmpi(out.modality_1,out.modality_2)
                                disp(['   The second ' out.modality_2 ' input will be matched to the ''aligned'' voxel space of the first for drawing mask VOI'])
                            else
                                disp(['   The ' out.modality_2 ' input will be matched to the ''aligned'' voxel space of the ' out.modality_1 ' for drawing mask VOI'])
                            end
                            disp(['   ...but output mask will be in voxel space of ' out.modality_1 ' input, i.e. not ''aligned'' space'])
                            struct1 = alignVolume(struct1);
                            disp('      aligning volume...')
                        else
                            if strcmpi(out.modality_1,out.modality_2)
                                disp(['   The second ' out.modality_2 ' input will be matched to the first ' out.modality_1 ' slice for drawing mask VOI'])
                            else
                                disp(['   The ' out.modality_2 ' input will be matched to the ' out.modality_1 ' slice for drawing mask VOI'])
                            end
                            if align
                                disp('   ...so ignoring ''align'' input')
                            end
                        end
                    end
                    temp = transformVoxelSpace(struct2,struct1);
                    if isempty(temp) || (isnumeric(temp) && temp==-1)
                        out = -1;
                        return
                    end
                    if isfield(struct1,'volume')
                        temp.vol_1 = struct1.volume;
                    elseif isfield(struct1,'volumes')
                        temp.vol_1 = struct1.volumes;
                    end
                    if isfield(temp,'volume')
                        temp.vol_2 = temp.volume;
                    elseif isfield(temp,'volumes')
                        temp.vol_2 = temp.volumes;
                    end
                else %extract 2D slice from 3d volume
                    %which is the 2D data?
                    if isfield(struct1,'volume')
                        siz1 = size(struct1.volume);
                    elseif isfield(struct1,'volumes')
                        siz1 = size(struct1.volumes);
                    end
                    if numel(siz1)==3
                        mod2D = out.modality_1;
                        mod3D = out.modality_2;
                    else
                        mod3D = out.modality_1;
                        mod2D = out.modality_2;
                    end                    
                    if strcmpi(mod2D,mod3D)
                        disp(['   The ' mod3D ' volume input will be matched to the other ' mod2D ' slice input'])
                    else
                        disp(['   The ' mod3D ' volume input will be matched to the ' mod2D ' slice input'])
                    end
                    if size(struct1.volumes,3)==1
                        temp = transformVoxelSpace(struct2,struct1);
                        if isempty(temp) || (isnumeric(temp) && temp==-1)
                            out = -1;
                            return
                        end
                        if isfield(struct1,'volume')
                            temp.vol_1 = struct1.volume;
                        elseif isfield(struct1,'volumes')
                            temp.vol_1 = struct1.volumes;
                        end
                        if isfield(temp,'volume')
                            temp.vol_2 = temp.volume;
                        elseif isfield(temp,'volumes')
                            temp.vol_2 = temp.volumes;
                        end
                    else
                        temp = transformVoxelSpace(struct1,struct2);
                        if isempty(temp) || (isnumeric(temp) && temp==-1)
                            out = -1;
                            return
                        end
                        if isfield(struct1,'volume')
                            temp.vol_2 = struct2.volume;
                        elseif isfield(struct1,'volumes')
                            temp.vol_2 = struct2.volumes;
                        end
                        if isfield(temp,'volume')
                            temp.vol_1 = temp.volume;
                        elseif isfield(temp,'volumes')
                            temp.vol_1 = temp.volumes;
                        end
                    end
                    if isfield(out,'getCoord') && out.getCoord
                        if numel(siz1)==3
                            if align
                                disp('   ...so ignoring ''align'' input')
                                disp(['   ...and output coordinates will be in voxel space of ' out.modality_1 ' volume, i.e. not slice-matched space'])
                            else
                                disp(['   ...but output coordinates will be in voxel space of ' out.modality_1 ' volume, i.e. not slice-matched space'])
                            end
                        end
                    end
                end
            end
            
            out.vol_1 = temp.vol_1;
            out.vol_2 = temp.vol_2;
            out.bedRange1 = single(temp.bedRange1);
            out.bedRange2 = single(temp.bedRange2);
            out.bedRange3 = single(temp.bedRange3);
            out.FOV1 = single(abs(diff(temp.bedRange1)));
            out.FOV2 = single(abs(diff(temp.bedRange2)));
            out.FOV3 = single(abs(diff(temp.bedRange3)));
            %get voxel absolute x, y, and z locations in mm
            out.locX = temp.locX;
            out.locY = temp.locY;
            out.locZ = temp.locZ;
            out.IOP = temp.IOP;
            
            if isfield(struct1,'units')
                out.units_1 = struct1.units;
                if SUVunits_1
                    out.vol_1 = out.vol_1 * SUVconv_1;
                    out.units_1 = 'SUV';
                end
            end
            if isfield(struct2,'units')
                out.units_2 = struct2.units;
                if SUVunits_2
                    out.vol_2 = out.vol_2 * SUVconv_2;
                    out.units_2 = 'SUV';
                end
            end
            if islogical(out.vol_1)% || all(out.vol_1(:)==0 | out.vol_1(:)==1)
                out.modality_1 = 'mask';
                if isfield(out,'units_1')
                    out = rmfield(out,'units_1');
                end
            end
            if islogical(out.vol_2)% || all(out.vol_2(:)==0 | out.vol_2(:)==1)
                out.modality_2 = 'mask';
                if isfield(out,'units_2')
                    out = rmfield(out,'units_2');
                end
            end
            if isfield(struct1,'patientOrientation')
                out.patientOrientation = struct1.patientOrientation;
            elseif isfield(struct2,'patientOrientation')
                out.patientOrientation = struct2.patientOrientation;
            end
        end
        out = checkMissingSlices(out);
        return
    end
end

end

function out = parseInputStruct(in)

fn = fieldnames(in);

for n=1:numel(fn)
    if strcmpi(fn{n},'vol') || strcmpi(fn{n},'vols') ||...
            strcmpi(fn{n},'volume') || strcmpi(fn{n},'volumes')
        if numel(size(in.(fn{n})))==4
            disp('Displaying last frame of dynamic series')
            %             disp(' ')
            out.vol = in.(fn{n})(:,:,:,end);
        else
            out.vol = in.(fn{n});
        end
    elseif strcmpi(fn{n},'modality')
        out.modality = in.(fn{n});
    elseif strcmpi(fn{n},'gantryModel')
        out.gantryModel = in.(fn{n});
    elseif strcmpi(fn{n},'FOV1') || strcmpi(fn{n},'FOV1')
        out.FOV1 = in.(fn{n});
    elseif strcmpi(fn{n},'bedRange1') || strcmpi(fn{n},'bedRangeX')
        out.bedRange1 = in.(fn{n});
        out.FOV1 = abs(diff(out.bedRange1));
    elseif strcmpi(fn{n},'FOV2') || strcmpi(fn{n},'FOV2')
        out.FOV2 = in.(fn{n});
    elseif strcmpi(fn{n},'bedRange2') || strcmpi(fn{n},'bedRangeY')
        out.bedRange2 = in.(fn{n});
        out.FOV2 = abs(diff(out.bedRange2));
    elseif strcmpi(fn{n},'FOV3') || strcmpi(fn{n},'FOV3')
        out.FOV3 = in.(fn{n});
    elseif strcmpi(fn{n},'bedRange3') || strcmpi(fn{n},'bedRangeZ')
        out.bedRange3 = in.(fn{n});
        out.FOV3 = abs(diff(out.bedRange3));
    elseif strcmpi(fn{n},'voxDim1')
        out.voxDim1 = in.voxDim1;
    elseif strcmpi(fn{n},'voxDim2')
        out.voxDim2 = in.voxDim2;
    elseif strcmpi(fn{n},'voxDim3')
        out.voxDim3 = in.voxDim3;
    elseif strcmpi(fn{n},'IOP')
        out.IOP = in.IOP;
    elseif strcmpi(fn{n},'IPP')
        out.IPP = in.IPP;
    elseif strcmpi(fn{n},'locX')
        out.locX = in.locX;
    elseif strcmpi(fn{n},'locY')
        out.locY = in.locY;
    elseif strcmpi(fn{n},'locZ')
        out.locZ = in.locZ;
    elseif strcmpi(fn{n},'filePath')
        out.filePath = in.(fn{n});
    elseif strcmpi(fn{n},'units')
        out.units = in.(fn{n});
    elseif strcmpi(fn{n},'patientOrientation')
        out.patientOrientation = in.(fn{n});
    end
end

end

function out = checkMissingSlices(in)

out = in;

if isstruct(in) && (isfield(in,'IPP') || isfield(in,'locX'))
    if ~isfield(in,'IPP')
        in.IPP = [squeeze(in.locX(1,1,1:end)) squeeze(in.locY(1,1,1:end)) squeeze(in.locZ(1,1,1:end))]';
    end
    %use slice positions
    temp = (in.IPP(:,2:end)-in.IPP(:,1:end-1))';
    znorms = sqrt(temp(:,1).^2+temp(:,2).^2+temp(:,3).^2);
    %find any gaps, this condition would fail if all missing slices were consistent regularly spaced
    idx = find(znorms>=1.9*min(znorms)); %allow 10% precision error
    if ~isempty(idx)
        disp('Filling in missing slices by interpolation')
        insertIdx = idx+1;
        nSlices2add = round(znorms(idx)/min(znorms))-1;        
        %make new IPP array
        newIPP = in.IPP(:,1:insertIdx(1)-1);
        for n=1:numel(insertIdx)
            %voxel dimensions
            xDim = (in.IPP(1,insertIdx(n))-in.IPP(1,insertIdx(n)-1))/(nSlices2add(n)+1);
            yDim = (in.IPP(2,insertIdx(n))-in.IPP(2,insertIdx(n)-1))/(nSlices2add(n)+1);
            zDim = (in.IPP(3,insertIdx(n))-in.IPP(3,insertIdx(n)-1))/(nSlices2add(n)+1);
            newIPP = [newIPP ...
                [in.IPP(1,insertIdx(n)-1)+xDim:xDim:in.IPP(1,insertIdx(n))-xDim;...
                in.IPP(2,insertIdx(n)-1)+yDim:yDim:in.IPP(2,insertIdx(n))-yDim;...
                in.IPP(3,insertIdx(n)-1)+zDim:zDim:in.IPP(3,insertIdx(n))-zDim]];
            if n==numel(insertIdx)
                newIPP = [newIPP in.IPP(:,insertIdx(n):end)];
            else
                newIPP = [newIPP in.IPP(:,insertIdx(n):insertIdx(n+1)-1)];
            end
        end
        %create template structures for transfom voxel space function
        temp1 = struct;
        temp2 = struct;
        
        if isfield(in,'vol')
            temp2.volumes = false(size(in.vol,1),size(in.vol,2),size(newIPP,2));
            temp2.imSz1 = size(in.vol,1);
            temp2.imSz2 = size(in.vol,2);
            temp2.imSz3 = size(newIPP,2);
            temp2.IOP = in.IOP;
            temp2.IPP = newIPP;
            temp2.voxDim1 = diff(in.bedRange1)/temp2.imSz1;
            temp2.voxDim2 = diff(in.bedRange2)/temp2.imSz2;
            temp2.voxDim3 = norm([xDim yDim zDim]);
            [temp2.bedRange1, temp2.bedRange2, temp2.bedRange3, ...
                temp2.locX, temp2.locY, temp2.locZ] = getVolumeBedRanges(temp2);
            temp1.volumes = in.vol;
            temp1.imSz1 = size(in.vol,1);
            temp1.imSz2 = size(in.vol,2);
            temp1.imSz3 = size(in.vol,3);
            temp1.IOP = in.IOP;
            temp1.locX = in.locX;
            temp1.locY = in.locY;
            temp1.locZ = in.locZ;
            
            temp = transformVoxelSpace(temp1,temp2);
            
            out.vol = temp.volumes;
            out.IPP = temp.IPP;
            out.bedRange1 = temp.bedRange1;
            out.bedRange2 = temp.bedRange2;
            out.bedRange3 = temp.bedRange3;
            out.FOV1 = single(abs(diff(temp.bedRange1)));
            out.FOV2 = single(abs(diff(temp.bedRange2)));
            out.FOV3 = single(abs(diff(temp.bedRange3)));
            out.locX = temp.locX;
            out.locY = temp.locY;
            out.locZ = temp.locZ;
        elseif isfield(in,'vol_1')
            temp2.volumes = false(size(in.vol_1,1),size(in.vol_1,2),size(newIPP,2));
            temp2.imSz1 = size(in.vol_1,1);
            temp2.imSz2 = size(in.vol_1,2);
            temp2.imSz3 = size(newIPP,2);
            temp2.IOP = in.IOP;
            temp2.IPP = newIPP;
            temp2.voxDim1 = diff(in.bedRange1)/temp2.imSz1;
            temp2.voxDim2 = diff(in.bedRange2)/temp2.imSz2;
            temp2.voxDim3 = norm([xDim yDim zDim]);
            [temp2.bedRange1, temp2.bedRange2, temp2.bedRange3, ...
                temp2.locX, temp2.locY, temp2.locZ] = getVolumeBedRanges(temp2);
            
            temp1.volumes = in.vol_1;
            temp1.imSz1 = size(in.vol_1,1);
            temp1.imSz2 = size(in.vol_1,2);
            temp1.imSz3 = size(in.vol_1,3);
            temp1.IOP = in.IOP;
            temp1.locX = in.locX;
            temp1.locY = in.locY;
            temp1.locZ = in.locZ;
            
            temp = transformVoxelSpace(temp1,temp2);
            
            out.vol_1 = temp.volumes;
            out.IPP = temp.IPP;
            out.bedRange1 = temp.bedRange1;
            out.bedRange2 = temp.bedRange2;
            out.bedRange3 = temp.bedRange3;
            out.FOV1 = single(abs(diff(temp.bedRange1)));
            out.FOV2 = single(abs(diff(temp.bedRange2)));
            out.FOV3 = single(abs(diff(temp.bedRange3)));
            out.locX = temp.locX;
            out.locY = temp.locY;
            out.locZ = temp.locZ;
            
            temp1.volumes = in.vol_2;
            temp = transformVoxelSpace(temp1,temp2);
            
            out.vol_2 = temp.volumes;        
        end
    end
end

end

function initialSetup(in)

global handles

if isfield(in,'FOV1')
    FOV1 = in.FOV1;
else
    FOV1 = in.voxDim1*size(in.vol,1);
end
if isfield(in,'FOV2')
    FOV2 = in.FOV2;
else
    FOV2 = in.voxDim2*size(in.vol,2);
end
if isfield(in,'FOV3')
    FOV3 = in.FOV3;
else
    FOV3 = in.voxDim3*size(in.vol,3);
end

screenPos = get(0,'ScreenSize');
screenSz = [screenPos(3) screenPos(4)];

%get scaling value for display text
handles.figure.DPIscl = 1;
s = settings;
try
    if s.matlab.desktop.HighDPISupport
        handles.figure.DPIscl = 1920/screenSz(1);
    end
catch
end
handles.figure.aspectRatio = 1.7778; %HD

if screenSz(1)/screenSz(2)>=handles.figure.aspectRatio
    figHeight = 0.8*screenSz(2);
    figSz = [figHeight*handles.figure.aspectRatio figHeight];
else
    figWidth = 0.8*screenSz(1);
    figSz = [figWidth figWidth/handles.figure.aspectRatio];
end
handles.figure.figPos = [(screenPos(3)-figSz(1))/2 (screenPos(4)-figSz(2))/2 figSz];
handles.figure.f1 = formatFigure(handles.figure.figPos);
makeInactive()

%TRANSAXIAL, CORONAL, AND SAGITTAL AXES POSITIONS
handles.figure.axGap = 0.035*figSz(1);

%DEFINE CENTER POINTS OF ALL AXES
handles.figure.ax1x = figSz(1)*0.2;
handles.figure.ax2x = figSz(1)*0.5;
handles.figure.ax3x = figSz(1)*0.8;
handles.figure.ax_y = figSz(2)*0.49;

%CALCULATE DISPLAY SCALE FACTOR AND AXES POSITIONS THEN FORMAT AXES
if FOV1>FOV3
    scl = (handles.figure.ax3x-handles.figure.ax1x-2*handles.figure.axGap)/(2*FOV1);
elseif FOV2>FOV3
    scl = (handles.figure.ax3x-handles.figure.ax1x-2*handles.figure.axGap)/(2*FOV2);
else
    scl = 0.5*figSz(2)/FOV3;
end
ax1Pos = [handles.figure.ax1x-scl*FOV1/2 handles.figure.ax_y-scl*FOV2/2 scl*FOV1 scl*FOV2];
ax2Pos = [handles.figure.ax2x-scl*FOV1/2 handles.figure.ax_y-scl*FOV3/2 scl*FOV1 scl*FOV3];
ax3Pos = [handles.figure.ax3x-scl*FOV2/2 handles.figure.ax_y-scl*FOV3/2 scl*FOV2 scl*FOV3];
while ax1Pos(1)+ax1Pos(3)>=ax2Pos(1)-handles.figure.axGap/2 || ...
        ax2Pos(1)+ax2Pos(3)>=ax3Pos(1)-handles.figure.axGap/2
    scl = scl-0.1;
    ax1Pos = [handles.figure.ax1x-scl*FOV1/2 handles.figure.ax_y-scl*FOV2/2 scl*FOV1 scl*FOV2];
    ax2Pos = [handles.figure.ax2x-scl*FOV1/2 handles.figure.ax_y-scl*FOV3/2 scl*FOV1 scl*FOV3];
    ax3Pos = [handles.figure.ax3x-scl*FOV2/2 handles.figure.ax_y-scl*FOV3/2 scl*FOV2 scl*FOV3];
end

[handles.axes.a1,handles.axes.a2,handles.axes.a3] = formatAxes(ax1Pos,ax2Pos,ax3Pos);
handles.axes.a1title = '';
handles.axes.a2title = '';
handles.axes.a3title = '';
set(handles.axes.a1,'ButtonDownFcn','')
set(handles.axes.a2,'ButtonDownFcn','')
set(handles.axes.a3,'ButtonDownFcn','')

%NOW SLIDER BARS FOR CONTROLLING SLICE NAVIGATION
handles.display.sliders.sldrWidth = 0.028*figSz(2);
%AXES SLIDER GAP
sldrGap = 0.025*figSz(1);
[handles.display.sliders.s1,handles.display.sliders.s2,handles.display.sliders.s3] = ...
    formatAxSliders(handles.display.sliders.sldrWidth,sldrGap,ax1Pos,ax2Pos,ax3Pos,FOV2,FOV3);

handles.figure.top_y = 0.88*figSz(2);

cleanupPreviousInstance()

if ~handles.display.images.fusedImages
    %IMAGE SMOOTHING DISPLAY PANEL
    smpSz = [0.22*figSz(1) 0.083*figSz(2)];
    smpPos = [handles.figure.ax1x-0.125*figSz(1) handles.figure.top_y-smpSz(2)/2 smpSz];
    s6Sz = [0.11*figSz(1) handles.display.sliders.sldrWidth];
    s6Pos = [smpPos(1)+(smpSz(1)-s6Sz(1))/2 smpPos(2)+handles.display.sliders.sldrWidth-0.006*figSz(1) s6Sz];
    [handles.display.sliders.s6,handles.display.sliders.s6t] = formatSmoothingPanel(smpPos,s6Pos,handles.figure.DPIscl);
    %COORDINATE DISPLAY PANEL
    cpSz = [0.105*figSz(1) 0.083*figSz(2)];
    % cpPos = [handles.figure.ax1x+0.12*figSz(1) handles.figure.top_y-cpSz(2)/2 cpSz];
    cpPos = [handles.figure.ax2x-cpSz(1)/2 handles.figure.top_y-cpSz(2)/2 cpSz];
    ctSz = [0.086*figSz(1) 0.042*figSz(2)];
    ctPos = [cpPos(1)+(cpSz(1)-ctSz(1))/2 cpPos(2)+0.004*figSz(2) ctSz];
    handles.display.text.cTxt = formatCoordPanel(cpPos,ctPos,handles.figure.DPIscl);
    %SLIDER BARS FOR CONTROLLING WINDOW DISPLAY RANGE
    spSz = [0.27*figSz(1) 2*handles.display.sliders.sldrWidth+0.07*figSz(2)];
    spPos = [handles.figure.ax3x-0.22*figSz(1) handles.figure.top_y-spSz(2)/2 spSz];
    sSz = [0.16*figSz(1) handles.display.sliders.sldrWidth];
    s4Pos = [spPos(1)+(spSz(1)-sSz(1))/2-0.015*figSz(1) spPos(2)+handles.display.sliders.sldrWidth+0.028*figSz(2) sSz];
    s5Pos = [spPos(1)+(spSz(1)-sSz(1))/2-0.015*figSz(1) spPos(2)+handles.display.sliders.sldrWidth-0.007*figSz(2) sSz];
    [handles.display.sliders.s4,handles.display.sliders.s5,handles.display.sliders.s4t,handles.display.sliders.s5t] = ...
        formatDispPanel(spPos,s4Pos,s5Pos,handles.figure.DPIscl);
    %COLORMAP (OR CT) PANEL
    cmSz = [0.055*figSz(1) 2*handles.display.sliders.sldrWidth+0.07*figSz(2)];
    cmPos = [handles.figure.ax3x+0.065*figSz(1) handles.figure.top_y-spSz(2)/2 cmSz];
    bSz = [(-1.4e-5*figSz(1)+0.055)*figSz(1)/handles.figure.DPIscl 0.03*figSz(2)];
    fs = 0.01 * figSz(1) * handles.figure.DPIscl;
    if handles.figure.DPIscl~=1
        fs = fs * 0.7;
    end
    if strcmpi(handles.volume.modality,'CT')
        formatCTPanel(cmPos,handles.figure.DPIscl)
        set(findobj('Tag','CMPanel'),'Title','Presets',...
            'FontSize',fs,'BorderType','etchedin')
    else
        formatCMPanel(cmPos,handles.figure.DPIscl)
        set(findobj('Tag','CMPanel'),'Title','Color',...
            'FontSize',fs,'BorderType','etchedin')
    end
    set(findobj('Tag','b1'),'Position',[cmPos(3)/2-(bSz(1)/2) 0.065*figSz(2) bSz])
    set(findobj('Tag','b2'),'Position',[cmPos(3)/2-(bSz(1)/2) 0.0375*figSz(2) bSz])
    set(findobj('Tag','b3'),'Position',[cmPos(3)/2-(bSz(1)/2) 0.01*figSz(2) bSz])
    %VOXEL VALUE TEXT BOX
    vpSz = [0.097*figSz(1) 0.083*figSz(2)];
    % vpPos = [handles.figure.ax2x-vpSz(1)/2 handles.figure.top_y-vpSz(2)/2 vpSz];
    vpPos = [handles.figure.ax1x+0.125*figSz(1) handles.figure.top_y-vpSz(2)/2 vpSz];
    vtSz = [0.094*figSz(1) 0.042*figSz(2)];
    vtPos = [vpPos(1)+(vpSz(1)-vtSz(1))/2 vpPos(2)+0.004*figSz(2) vtSz];
    handles.display.text.vTxt = formatVoxPanel(vpPos,vtPos,handles.figure.DPIscl);
    %use quanitfication units if defined
    if ~strcmpi(handles.display.images.units,'')
        switch lower(handles.display.images.units)
            case 'bqml'
                actTxt = 'Bq/mL';
            case 'hu'
                actTxt = 'HU';
            case 'suv'
                actTxt = 'SUV';
            case 'propcps'
                actTxt = 'CPS';
            case 'cnts'
                actTxt = 'CNTS';
            case '1/cm'
                actTxt = '1/cm';
            case 'ml/ccm/min'
                actTxt = 'mL/ccm/min';
            otherwise
                warning(['String is not defined for activity units ''' handles.display.images.units ''''])
                actTxt = '(undef)';
        end
        set(findobj('Tag','VoxPanel'),'Title',actTxt)
    end
else
    %image weights for fusion
    fwpSz = [0.23*figSz(1) 0.083*figSz(2)];
    fwpPos = [handles.figure.ax1x-0.125*figSz(1) handles.figure.top_y-fwpSz(2)/2 fwpSz];
    s6Sz = [0.15*figSz(1) handles.display.sliders.sldrWidth];
    s6Pos = [fwpPos(1)+(fwpSz(1)-s6Sz(1))/2 fwpPos(2)+handles.display.sliders.sldrWidth-0.006*figSz(1) s6Sz];
    [handles.display.sliders.s6_1] = formatFusedWeightPanel(fwpPos,s6Pos,handles.figure.DPIscl);
    set(findobj('Tag','s6Text_1'),'String',handles.volume.modality_1)
    set(findobj('Tag','s6Text_2'),'String',handles.volume.modality_2)
    %slider bars for controlling display ranges for each image
    spSz = [0.47*figSz(1) 2*handles.display.sliders.sldrWidth+0.08*figSz(2)];
    spPos = [handles.figure.ax3x-0.35*figSz(1) handles.figure.top_y-spSz(2)/2 spSz];
    %for 2 panels (both imaages)
    sSz = [0.09*figSz(1) handles.display.sliders.sldrWidth];
    s4Pos1 = [spPos(1)+0.040*figSz(1) spPos(2)+handles.display.sliders.sldrWidth+0.028*figSz(2) sSz];
    s4Pos2 = [spPos(1)+0.275*figSz(1) spPos(2)+handles.display.sliders.sldrWidth+0.028*figSz(2) sSz];
    s5Pos1 = [spPos(1)+0.040*figSz(1) spPos(2)+handles.display.sliders.sldrWidth-0.007*figSz(2) sSz];
    s5Pos2 = [spPos(1)+0.275*figSz(1) spPos(2)+handles.display.sliders.sldrWidth-0.007*figSz(2) sSz];
    [handles.display.sliders.s4_1,handles.display.sliders.s4_2,handles.display.sliders.s5_1,handles.display.sliders.s5_2,...
        handles.display.sliders.s4t_1,handles.display.sliders.s4t_2,handles.display.sliders.s5t_1,handles.display.sliders.s5t_2] = ...
        formatFusedDispPanel(spPos,s4Pos1,s4Pos2,s5Pos1,s5Pos2,handles.figure.DPIscl);
    set(findobj('Tag','imLabel_1'),'String',handles.volume.modality_1)
    set(findobj('Tag','imLabel_2'),'String',handles.volume.modality_2)
    %CM OR CT PRESET PANELS
    cmSz = [0.035*figSz(1) 2*handles.display.sliders.sldrWidth+0.045*figSz(2)];
    cmPos_1 = [s4Pos1(1)+s4Pos1(3)+0.067*figSz(1) handles.figure.top_y-spSz(2)/2.15 cmSz];
    cmPos_2 = [s4Pos2(1)+s4Pos2(3)+0.067*figSz(1) handles.figure.top_y-spSz(2)/2.15 cmSz];
    if strcmpi(handles.volume.modality_1,'CT')
        formatCTPanel(cmPos_1,handles.figure.DPIscl,...
            'CMPanel_1','b1_1','b2_1','b3_1')
    else
        formatCMPanel(cmPos_1,handles.figure.DPIscl,...
            'CMPanel_1','b1_1','b2_1','b3_1')
    end
    if strcmpi(handles.volume.modality_2,'CT')
        formatCTPanel(cmPos_2,handles.figure.DPIscl,...
            'CMPanel_2','b1_2','b2_2','b3_2')
    else
        formatCMPanel(cmPos_2,handles.figure.DPIscl,...
            'CMPanel_2','b1_2','b2_2','b3_2')
    end
    %voxel values text box for PET and CT
    vpSz = [0.12*figSz(1) 2*handles.display.sliders.sldrWidth+0.08*figSz(2)];
    vpPos = [handles.figure.ax1x+0.117*figSz(1) handles.figure.top_y-vpSz(2)/2 vpSz];
    vtSz = [0.07*figSz(1) handles.display.sliders.sldrWidth];
    vtPos1 = [vpPos(1)+(vpSz(1)-vtSz(1))*0.85 vpPos(2)+handles.display.sliders.sldrWidth+0.028*figSz(2) vtSz];
    vtPos2 = [vpPos(1)+(vpSz(1)-vtSz(1))*0.85 vpPos(2)+handles.display.sliders.sldrWidth-0.007*figSz(2) vtSz];
    [handles.display.text.vTxt_1,handles.display.text.vTxt_2] = ...
        formatFusedVoxPanel(vpPos,vtPos1,vtPos2,handles.figure.DPIscl);
    set(findobj('Tag','imVoxLabel_1'),'String',handles.volume.modality_1)
    set(findobj('Tag','imVoxLabel_2'),'String',handles.volume.modality_2)
    %use quanitfication units if defined
    if ~strcmpi(handles.display.images.units_1,'')
        switch lower(handles.display.images.units_1)
            case 'bqml'
                actTxt = 'Bq/mL';
            case 'hu'
                actTxt = 'HU';
            case 'suv'
                actTxt = 'SUV';
            case 'propcps'
                actTxt = 'CPS';
            case 'cnts'
                actTxt = 'CNTS';
            case '1/cm'
                actTxt = '1/cm';
            case 'ml/ccm/min'
                actTxt = 'mL/ccm/min';
            otherwise
                warning(['String is not defined for activity units ''' handles.display.images.units_1 ''''])
                actTxt = '(undef)';
        end
        set(findobj('Tag','imVoxLabel_1'),'String',actTxt)
    end
    if ~strcmpi(handles.display.images.units_2,'')
        switch lower(handles.display.images.units_2)
            case 'bqml'
                actTxt = 'Bq/mL';
            case 'hu'
                actTxt = 'HU';
            case 'suv'
                actTxt = 'SUV';
            case 'propcps'
                actTxt = 'CPS';
            case 'cnts'
                actTxt = 'CNTS';
            case '1/cm'
                actTxt = '1/cm';
            case 'ml/ccm/min'
                actTxt = 'mL/ccm/min';
            otherwise
                warning(['String is not defined for activity units ''' handles.display.images.units_2 ''''])
                actTxt = '(undef)';
        end
        set(findobj('Tag','imVoxLabel_2'),'String',actTxt)
    end
end

handles.display.sliders.sliderGrab = false;

txtStr = {'{''Scroll to zoom on current position''}'};
%java
labelStr = {['<html><div style="text-align: center;"><font size="' num2str(0.007*figSz(1)*handles.figure.DPIscl) ...
    '"><b>Scroll</b> to zoom on current position<br></font>']};
txtCnt = 1;
if ~handles.display.images.fusedImages && ~strcmpi(handles.volume.modality,'CT')
    txtStr = {[txtStr{:} ' {''Double click to invert display colormap''}']};
    labelStr = {[labelStr{:} ['<font size="' num2str(0.007*figSz(1)*handles.figure.DPIscl) ...
        '"><b>Double click</b> to invert display colormap<br></font>']]};
    txtCnt = txtCnt+1;
end
if handles.display.images.drawROI
    txtStr = [txtStr {' {''Right click to fill in enclosed boundary or remove pixel blocks''}'}];
    labelStr = [labelStr {'<font size="' num2str(0.007*figSz(1)*handles.figure.DPIscl) ...
        '"><b>Right click</b> to fill in enclosed boundary or remove pixel blocks<br></font>'}];
    txtCnt = txtCnt+1;
end

ctSz = [0.4*figSz(1) txtCnt*0.025*figSz(2)];
ctPos = [(figSz(1)-ctSz(1))/2 0.02*figSz(2) ctSz];

%use Java swing control if possible for text customization
try
    import javax.swing.JLabel;
    jLabel = javaObjectEDT('javax.swing.JLabel',[labelStr{1} '</html>']);
    [handles.display.text.instructions,handles.display.text.instructionsWrapper] = ...
        javacomponent(jLabel,ctPos,gcf);
    handles.display.text.instructions.setHorizontalAlignment(JLabel.CENTER);
    handles.display.text.instructions.setVerticalAlignment(JLabel.TOP);
    handles.display.text.instructionsStr = labelStr;
    handles.display.text.instructionsJavaBool = true;
    set(handles.display.text.instructionsWrapper,'UserData',handles.display.text.instructions)
    set(handles.display.text.instructionsWrapper,'Tag','InstructionsText')
catch %or use default matlab uicontrol
    str = ['''Style'',''Text'','...
        '''Units'',''Pixel'','...
        '''Position'',[' num2str(ctPos) '],'...
        '''FontSize'',' num2str(0.007*figSz(1)*handles.figure.DPIscl) ','...
        '''FontWeight'',''Normal'','...
        '''String'',[' txtStr{1} '],'...
        '''BackgroundColor'',[' num2str(get(findobj('Tag','3DviewerFig'),'color')) ']'];
    handles.display.text.instructions = initiateObject('InstructionsText','uicontrol',str);
    uistack(handles.display.text.instructions,'bottom')
    handles.display.text.instructionsStr = [eval(['[' txtStr{:} ']'])];
    handles.display.text.instructionsJavaBool = false;
end

end

function cleanupPreviousInstance

%from single image mode
delete(findobj('Tag','VoxText'))
delete(findobj('Tag','SmoothingPanel'))
delete(findobj('Tag','s6Slider'))
delete(findobj('Tag','s6Text1'))
delete(findobj('Tag','s6Text2'))
delete(findobj('Tag','s6Text2units'))
delete(findobj('Tag','s6Text_1'))
delete(findobj('Tag','s6Text_2'))
delete(findobj('Tag','CoordinatePanel'))
delete(findobj('Tag','CoordinateText'))
delete(findobj('Tag','s4Slider'))
delete(findobj('Tag','s4Button1'))
delete(findobj('Tag','s4Text2'))
delete(findobj('Tag','s5Slider'))
delete(findobj('Tag','s5Button1'))
delete(findobj('Tag','s5Text2'))
delete(findobj('Tag','CMPanel'))
delete(findobj('Tag','b1'))
delete(findobj('Tag','b2'))
delete(findobj('Tag','b3'))
%from fused image mode
delete(findobj('Tag','FusedWeightPanel'))
delete(findobj('Tag','VoxText_1'))
delete(findobj('Tag','VoxText_2'))
delete(findobj('Tag','imVoxLabel_1'))
delete(findobj('Tag','imVoxLabel_2'))
delete(findobj('Tag','s6Slider'))
delete(findobj('Tag','s6Text1'))
delete(findobj('Tag','s6Text2'))
delete(findobj('Tag','s4Slider_1'))
delete(findobj('Tag','s4Slider_2'))
delete(findobj('Tag','s5Slider_1'))
delete(findobj('Tag','s5Slider_2'))
delete(findobj('Tag','s4Button1_1'))
delete(findobj('Tag','s4Text2_1'))
delete(findobj('Tag','s5Button1_1'))
delete(findobj('Tag','s5Text2_1'))
delete(findobj('Tag','s4Button1_2'))
delete(findobj('Tag','s4Text2_2'))
delete(findobj('Tag','s5Button1_2'))
delete(findobj('Tag','s5Text2_2'))
delete(findobj('Tag','CMPanel_1'))
delete(findobj('Tag','CMPanel_2'))
delete(findobj('Tag','imLabel_1'))
delete(findobj('Tag','imLabel_2'))
delete(findobj('Tag','b1_1'))
delete(findobj('Tag','b2_1'))
delete(findobj('Tag','b3_1'))
delete(findobj('Tag','b1_2'))
delete(findobj('Tag','b2_2'))
delete(findobj('Tag','b3_2'))
%always delete crosshairs so they are initiated on top
delete(findobj('Tag','a1crossH1'))
delete(findobj('Tag','a1crossH2'))
delete(findobj('Tag','a1crossV1'))
delete(findobj('Tag','a1crossV2'))
delete(findobj('Tag','a2crossH1'))
delete(findobj('Tag','a2crossH2'))
delete(findobj('Tag','a2crossV1'))
delete(findobj('Tag','a2crossV2'))
delete(findobj('Tag','a3crossH1'))
delete(findobj('Tag','a3crossH2'))
delete(findobj('Tag','a3crossV1'))
delete(findobj('Tag','a3crossV2'))

delete(findobj('Tag','LocText1x'))
delete(findobj('Tag','LocText1y'))
delete(findobj('Tag','LocText1z'))
delete(findobj('Tag','LocText2x'))
delete(findobj('Tag','LocText2y'))
delete(findobj('Tag','LocText2z'))
%just in case
delete(findobj('Tag','drawBtn'))
delete(findobj('Tag','doneBtn'))
delete(findobj('Tag','maskPatch1'))
delete(findobj('Tag','maskPatch2'))
delete(findobj('Tag','maskPatch3'))

try
    %for Java swing component
    temp = findobj('Tag','InstructionsText');
    delete(get(temp,'UserData'))
catch
    delete(findobj('Tag','InstructionsText'))
end

end

function f = formatFigure(figPos)

test = findobj('Tag','3DviewerFig');
if ~isempty(test)
    figure(test)
    pos = get(test,'Position');
    cntrPnt = [pos(1)+pos(3)/2 pos(2)+pos(4)/2];
    figPos = [cntrPnt(1)-figPos(3)/2 cntrPnt(2)-figPos(4)/2 figPos(3) figPos(4)];
end

str = ['''Name'',''3D Viewer'','...
    '''Units'',''Pixel'','...
    '''Position'',[' num2str(figPos) '],'...
    '''MenuBar'',''None'','...
    '''ToolBar'',''None'','...
    '''Resize'',''Off'','...
    '''NumberTitle'',''Off'''];
f = initiateObject('3DviewerFig','figure',str);

end

function [ax1,ax2,ax3] = formatAxes(ax1pos,ax2pos,ax3pos)


str = ['''Units'',''Pixel'','...
    '''Position'',[' num2str(ax1pos) ']'];
ax1 = initiateObject('axes1','axes',str);
box on; hold on

str = ['''Units'',''Pixel'','...
    '''Position'',[' num2str(ax2pos) ']'];
ax2 = initiateObject('axes2','axes',str);
box on; hold on

str = ['''Units'',''Pixel'','...
    '''Position'',[' num2str(ax3pos) ']'];
ax3 = initiateObject('axes3','axes',str);
box on; hold on

set(ax1,'XTick',[],'YTick',[])
set(ax2,'XTick',[],'YTick',[])
set(ax3,'XTick',[],'YTick',[])

end

function [s1,s2,s3] = formatAxSliders(sldrWidth,sldrGap,ax1Pos,ax2Pos,ax3Pos,FOV2,FOV3)

if FOV2>FOV3 %FOV1 IS NOT REPRESENTED VERTICALLY
    s1Pos = [ax1Pos(1) ax1Pos(2)-sldrWidth-sldrGap ax1Pos(3) sldrWidth];
    s2Pos = [ax2Pos(1) ax1Pos(2)-sldrWidth-sldrGap ax2Pos(3) sldrWidth];
    s3Pos = [ax3Pos(1) ax1Pos(2)-sldrWidth-sldrGap ax3Pos(3) sldrWidth];
else
    s1Pos = [ax1Pos(1) ax2Pos(2)-sldrWidth-sldrGap ax1Pos(3) sldrWidth];
    s2Pos = [ax2Pos(1) ax2Pos(2)-sldrWidth-sldrGap ax2Pos(3) sldrWidth];
    s3Pos = [ax3Pos(1) ax2Pos(2)-sldrWidth-sldrGap ax3Pos(3) sldrWidth];
end
%IF SLIDERS FOUND, DELETE TO REMOVE LISTENERS
delete(findobj('Tag','s1Slider'))
delete(findobj('Tag','s2Slider'))
delete(findobj('Tag','s3Slider'))
s1 = uicontrol('Style','Slider','Position',s1Pos,'Enable','off','Tag','s1Slider');
s2 = uicontrol('Style','Slider','Position',s2Pos,'Enable','off','Tag','s2Slider');
s3 = uicontrol('Style','Slider','Position',s3Pos,'Enable','off','Tag','s3Slider');

end

function vt = formatVoxPanel(vpPos,vtPos,txtScl)

if txtScl~=1 %high DPI scaling
    txtScl = txtScl * 0.7;
end

figSz = get(findobj('Tag','3DviewerFig'),'Position');
figSz = [figSz(3) figSz(4)];

str = ['''Title'',''Voxel Value'','...
    '''Units'',''Pixel'','...
    '''Position'',[' num2str(vpPos) '],'...
    '''FontSize'',' num2str(0.01*figSz(1)*txtScl) ','...
    '''TitlePosition'',''CenterTop'','...
    '''BackgroundColor'',[' num2str(get(findobj('Tag','3DviewerFig'),'color')) ']'];
initiateObject('VoxPanel','uipanel',str);

str = ['''Style'',''Text'','...
    '''Units'',''Pixel'','...
    '''Position'',[' num2str(vtPos) '],'...
    '''FontSize'',' num2str(0.01*figSz(1)*txtScl) ','...
    '''FontWeight'',''Bold'','...
    '''HorizontalAlignment'',''Center'','...
    '''BackgroundColor'',[' num2str(get(findobj('Tag','3DviewerFig'),'color')) ']'];
vt = initiateObject('VoxText','uicontrol',str);

end

function [vt1,vt2] = formatFusedVoxPanel(vpPos,vtPos1,vtPos2,txtScl)

if txtScl~=1 %high DPI scaling
    txtScl = txtScl * 0.7;
end

figSz = get(findobj('Tag','3DviewerFig'),'Position');
figSz = [figSz(3) figSz(4)];

str = ['''Title'',''Voxel Value'','...
    '''Units'',''Pixel'','...
    '''Position'',[' num2str(vpPos) '],'...
    '''FontSize'',' num2str(0.01*figSz(1)*txtScl) ','...
    '''TitlePosition'',''CenterTop'','...
    '''BackgroundColor'',[' num2str(get(findobj('Tag','3DviewerFig'),'color')) ']'];
initiateObject('VoxPanel','uipanel',str);

str = ['''Style'',''Text'','...
    '''Units'',''Pixel'','...
    '''Position'',[' num2str(vtPos1) '],'...
    '''FontSize'',' num2str(0.01*figSz(1)*txtScl) ','...
    '''FontWeight'',''Bold'','...
    '''HorizontalAlignment'',''Center'','...
    '''BackgroundColor'',[' num2str(get(findobj('Tag','3DviewerFig'),'color')) ']'];
vt1 = initiateObject('VoxText_1','uicontrol',str);

str = ['''Style'',''Text'','...
    '''Units'',''Pixel'','...
    '''Position'',[' num2str(vtPos2) '],'...
    '''FontSize'',' num2str(0.01*figSz(1)*txtScl) ','...
    '''FontWeight'',''Bold'','...
    '''HorizontalAlignment'',''Center'','...
    '''BackgroundColor'',[' num2str(get(findobj('Tag','3DviewerFig'),'color')) ']'];
vt2 = initiateObject('VoxText_2','uicontrol',str);

%add labels for each voxel value
txtWdth = 0.04*figSz(1);
txtHght = vtPos1(4);
pos1 = [vpPos(1)+(0.012*figSz(1)) vtPos1(2) txtWdth txtHght];
pos2 = [vpPos(1)+(0.012*figSz(1)) vtPos2(2) txtWdth txtHght];
im1str = '';
im2str = '';

str = ['''Style'',''Text'','...
    '''Units'',''Pixel'','...
    '''Position'',[' num2str(pos1) '],'...
    '''FontSize'',' num2str(0.01*figSz(1)*txtScl) ','...
    '''FontWeight'',''Normal'','...
    '''HorizontalAlignment'',''Center'','...
    '''String'',''' im1str ''','...
    '''BackgroundColor'',[' num2str(get(findobj('Tag','3DviewerFig'),'color')) ']'];
initiateObject('imVoxLabel_1','uicontrol',str);

str = ['''Style'',''Text'','...
    '''Units'',''Pixel'','...
    '''Position'',[' num2str(pos2) '],'...
    '''FontSize'',' num2str(0.01*figSz(1)*txtScl) ','...
    '''FontWeight'',''Normal'','...
    '''HorizontalAlignment'',''Center'','...
    '''String'',''' im2str ''','...
    '''BackgroundColor'',[' num2str(get(findobj('Tag','3DviewerFig'),'color')) ']'];
initiateObject('imVoxLabel_2','uicontrol',str);

end

function ct = formatCoordPanel(cpPos,ctPos,txtScl)

if txtScl~=1 %high DPI scaling
    txtScl = txtScl * 0.7;
end

figSz = get(findobj('Tag','3DviewerFig'),'Position');
figSz = [figSz(3) figSz(4)];

str = ['''Title'',''Coordinates'','...
    '''Units'',''Pixel'','...
    '''Position'',[' num2str(cpPos) '],'...
    '''FontSize'',' num2str(0.01*figSz(1)*txtScl) ','...
    '''TitlePosition'',''CenterTop'','...
    '''BackgroundColor'',[' num2str(get(findobj('Tag','3DviewerFig'),'color')) ']'];
initiateObject('CoordinatePanel','uipanel',str);

str = ['''Style'',''Text'','...
    '''Units'',''Pixel'','...
    '''Position'',[' num2str(ctPos) '],'...
    '''FontSize'',' num2str(0.01*figSz(1)*txtScl) ','...
    '''FontWeight'',''Bold'','...
    '''HorizontalAlignment'',''Center'','...
    '''BackgroundColor'',[' num2str(get(findobj('Tag','3DviewerFig'),'color')) ']'];
ct = initiateObject('CoordinateText','uicontrol',str);

end

function [lt2x,lt2y,lt2z] = formatLocPanel

global handles

figSz = get(handles.figure.f1,'Position');
figSz = [figSz(3) figSz(4)];

txtStr1x = '''x:''';
txtStr1y = '''y:''';
txtStr1z = '''z:''';
txtStr2x = '''''';
txtStr2y = '''''';
txtStr2z = '''''';

fs = 0.007*figSz(1)*handles.figure.DPIscl;

t1Sz = [0.02*figSz(1) 0.03*figSz(2)];
t1xPos = [0.005*figSz(1) 0.008*figSz(2)*fs t1Sz];
t1yPos = [0.005*figSz(1) 0.0055*figSz(2)*fs t1Sz];
t1zPos = [0.005*figSz(1) 0.003*figSz(2)*fs t1Sz];
t2Sz = [0.07*figSz(1) 0.03*figSz(2)];
t2xPos = [t1xPos(1)+t1Sz(1)+0.005*figSz(1) 0.008*figSz(2)*fs t2Sz];
t2yPos = [t1yPos(1)+t1Sz(1)+0.005*figSz(1) 0.0055*figSz(2)*fs t2Sz];
t2zPos = [t1zPos(1)+t1Sz(1)+0.005*figSz(1) 0.003*figSz(2)*fs t2Sz];

str = ['''Style'',''Text'','...
    '''Units'',''Pixel'','...
    '''Position'',[' num2str(t1xPos) '],'...
    '''FontSize'',' num2str(fs) ','...
    '''HorizontalAlignment'',''Right'','...
    '''FontWeight'',''Normal'','...
    '''String'',' txtStr1x ','...
    '''BackgroundColor'',[' num2str(get(findobj('Tag','3DviewerFig'),'color')) ']'];
lt1x = initiateObject('LocText1x','uicontrol',str);
str = ['''Style'',''Text'','...
    '''Units'',''Pixel'','...
    '''Position'',[' num2str(t1yPos) '],'...
    '''FontSize'',' num2str(fs) ','...
    '''HorizontalAlignment'',''Right'','...
    '''FontWeight'',''Normal'','...
    '''String'',' txtStr1y ','...
    '''BackgroundColor'',[' num2str(get(findobj('Tag','3DviewerFig'),'color')) ']'];
lt1y = initiateObject('LocText1y','uicontrol',str);
str = ['''Style'',''Text'','...
    '''Units'',''Pixel'','...
    '''Position'',[' num2str(t1zPos) '],'...
    '''FontSize'',' num2str(fs) ','...
    '''HorizontalAlignment'',''Right'','...
    '''FontWeight'',''Normal'','...
    '''String'',' txtStr1z ','...
    '''BackgroundColor'',[' num2str(get(findobj('Tag','3DviewerFig'),'color')) ']'];
lt1z = initiateObject('LocText1z','uicontrol',str);

str = ['''Style'',''Text'','...
    '''Units'',''Pixel'','...
    '''Position'',[' num2str(t2xPos) '],'...
    '''FontSize'',' num2str(fs) ','...
    '''HorizontalAlignment'',''Left'','...
    '''FontWeight'',''Normal'','...
    '''String'',' txtStr2x ','...
    '''BackgroundColor'',[' num2str(get(findobj('Tag','3DviewerFig'),'color')) ']'];
lt2x = initiateObject('LocText2x','uicontrol',str);
str = ['''Style'',''Text'','...
    '''Units'',''Pixel'','...
    '''Position'',[' num2str(t2yPos) '],'...
    '''FontSize'',' num2str(fs) ','...
    '''HorizontalAlignment'',''Left'','...
    '''FontWeight'',''Normal'','...
    '''String'',' txtStr2y ','...
    '''BackgroundColor'',[' num2str(get(findobj('Tag','3DviewerFig'),'color')) ']'];
lt2y = initiateObject('LocText2y','uicontrol',str);
str = ['''Style'',''Text'','...
    '''Units'',''Pixel'','...
    '''Position'',[' num2str(t2zPos) '],'...
    '''FontSize'',' num2str(fs) ','...
    '''HorizontalAlignment'',''Left'','...
    '''FontWeight'',''Normal'','...
    '''String'',' txtStr2z ','...
    '''BackgroundColor'',[' num2str(get(findobj('Tag','3DviewerFig'),'color')) ']'];
lt2z = initiateObject('LocText2z','uicontrol',str);

end

function [s4,s5,s4t,s5t] = formatDispPanel(spPos,s4Pos,s5Pos,txtScl)

if txtScl~=1 %high DPI scaling
    txtScl = txtScl * 0.7;
end

figSz = get(findobj('Tag','3DviewerFig'),'Position');
figSz = [figSz(3) figSz(4)];

%SLIDER PANEL
str = ['''Title'',''Viewing Window'','...
    '''Units'',''Pixel'','...
    '''Position'',[' num2str(spPos) '],'...
    '''FontSize'',' num2str(0.01*figSz(1)*txtScl) ','...
    '''TitlePosition'',''CenterTop'','...
    '''BackgroundColor'',[' num2str(get(findobj('Tag','3DviewerFig'),'color')) ']'];
initiateObject('DisplayPanel','uipanel',str);
%S4 SLIDER
delete(findobj('Tag','s4Slider'));
s4 = uicontrol('Style','Slider','Position',s4Pos,'Enable','off','Tag','s4Slider');
%S4 TEXT 1
str = ['''Style'',''PushButton'','...
    '''Units'',''Pixel'','...
    '''Position'',[' num2str([spPos(1)+0.006*figSz(1) s4Pos(2) 0.03*figSz(1) s4Pos(4)]) '],'...
    '''FontSize'',' num2str(0.01*figSz(1)*txtScl) ','...
    '''HorizontalAlignment'',''Center'','...
    '''String'',''Max'','...
    '''Callback'',@maxDispButtonCallback,'...
    '''BackgroundColor'',[' num2str(get(findobj('Tag','3DviewerFig'),'color')) ']'];
initiateObject('s4Button1','uicontrol',str);
%S4 TEXT 2
str = ['''Style'',''Edit'','...
    '''Units'',''Pixel'','...
    '''Position'',[' num2str([s4Pos(1)+s4Pos(3)+0.0092*figSz(1) s4Pos(2) 0.05*figSz(1) s4Pos(4)]) '],'...
    '''FontSize'',' num2str(0.01*figSz(1)*txtScl) ','...
    '''FontWeight'',''Bold'','...
    '''HorizontalAlignment'',''Center'','...
    '''String'','''','...
    '''Callback'',@maxDispTextCallback,'...
    '''BackgroundColor'',[' num2str(get(findobj('Tag','3DviewerFig'),'color')) ']'];
s4t = initiateObject('s4Text2','uicontrol',str);
%S5 SLIDER
delete(findobj('Tag','s5Slider'))
s5 = uicontrol('Style','Slider','Position',s5Pos,'Enable','off','Tag','s5Slider');
%S5 TEXT 1
str = ['''Style'',''PushButton'','...
    '''Units'',''Pixel'','...
    '''Position'',[' num2str([spPos(1)+0.006*figSz(1) s5Pos(2) 0.03*figSz(1) s5Pos(4)]) '],'...
    '''FontSize'',' num2str(0.01*figSz(1)*txtScl) ','...
    '''HorizontalAlignment'',''Center'','...
    '''String'',''Min'','...
    '''Callback'',@minDispButtonCallback,'...
    '''BackgroundColor'',[' num2str(get(findobj('Tag','3DviewerFig'),'color')) ']'];
initiateObject('s5Button1','uicontrol',str);
%S5 TEXT 2
str = ['''Style'',''Edit'','...
    '''Units'',''Pixel'','...
    '''Position'',[' num2str([s5Pos(1)+s5Pos(3)+0.0092*figSz(1) s5Pos(2) 0.05*figSz(1) s5Pos(4)]) '],'...
    '''FontSize'',' num2str(0.01*figSz(1)*txtScl) ','...
    '''FontWeight'',''Bold'','...
    '''HorizontalAlignment'',''Center'','...
    '''String'','''','...
    '''Callback'',@minDispTextCallback,'...
    '''BackgroundColor'',[' num2str(get(findobj('Tag','3DviewerFig'),'color')) ']'];
s5t = initiateObject('s5Text2','uicontrol',str);

end

function [s41,s42,s51,s52,s4t1,s4t2,s5t1,s5t2] = formatFusedDispPanel(spPos,s4Pos1,s4Pos2,s5Pos1,s5Pos2,txtScl)

if txtScl~=1 %high DPI scaling
    txtScl = txtScl * 0.7;
end

figSz = get(findobj('Tag','3DviewerFig'),'Position');
figSz = [figSz(3) figSz(4)];

%SLIDER PANEL
str = ['''Title'',''Fused Viewing Window'','...
    '''Units'',''Pixel'','...
    '''Position'',[' num2str(spPos) '],'...
    '''FontSize'',' num2str(0.01*figSz(1)*txtScl) ','...
    '''TitlePosition'',''CenterTop'','...
    '''BackgroundColor'',[' num2str(get(findobj('Tag','3DviewerFig'),'color')) ']'];
initiateObject('DisplayPanel','uipanel',str);
%S41 SLIDER
delete(findobj('Tag','s4Slider_1'));
s41 = uicontrol('Style','Slider','Position',s4Pos1,'Enable','off','Tag','s4Slider_1');
%S41 TEXT 1
str = ['''Style'',''PushButton'','...
    '''Units'',''Pixel'','...
    '''Position'',[' num2str([s4Pos1(1)-0.034*figSz(1) s4Pos1(2) 0.03*figSz(1) s4Pos1(4)]) '],'...
    '''FontSize'',' num2str(0.01*figSz(1)*txtScl) ','...
    '''HorizontalAlignment'',''Center'','...
    '''String'',''Max'','...
    '''Callback'',@maxDispButtonCallback_1,'...
    '''BackgroundColor'',[' num2str(get(findobj('Tag','3DviewerFig'),'color')) ']'];
initiateObject('s4Button1_1','uicontrol',str);
%S41 TEXT 2
str = ['''Style'',''Edit'','...
    '''Units'',''Pixel'','...
    '''Position'',[' num2str([s4Pos1(1)+s4Pos1(3)+0.008*figSz(1) s4Pos1(2) 0.05*figSz(1) s4Pos1(4)]) '],'...
    '''FontSize'',' num2str(0.01*figSz(1)*txtScl) ','...
    '''FontWeight'',''Bold'','...
    '''HorizontalAlignment'',''Center'','...
    '''String'','''','...
    '''Callback'',@maxDispTextCallback_1,'...
    '''BackgroundColor'',[' num2str(get(findobj('Tag','3DviewerFig'),'color')) ']'];
s4t1 = initiateObject('s4Text2_1','uicontrol',str);
%S42 SLIDER
delete(findobj('Tag','s4Slider_2'));
s42 = uicontrol('Style','Slider','Position',s4Pos2,'Enable','off','Tag','s4Slider_2');
%S42 TEXT 1
str = ['''Style'',''PushButton'','...
    '''Units'',''Pixel'','...
    '''Position'',[' num2str([s4Pos2(1)-0.034*figSz(1) s4Pos2(2) 0.03*figSz(1) s4Pos2(4)]) '],'...
    '''FontSize'',' num2str(0.01*figSz(1)*txtScl) ','...
    '''HorizontalAlignment'',''Center'','...
    '''String'',''Max'','...
    '''Callback'',@maxDispButtonCallback_2,'...
    '''BackgroundColor'',[' num2str(get(findobj('Tag','3DviewerFig'),'color')) ']'];
initiateObject('s4Button1_2','uicontrol',str);
%S42 TEXT 2
str = ['''Style'',''Edit'','...
    '''Units'',''Pixel'','...
    '''Position'',[' num2str([s4Pos2(1)+s4Pos2(3)+0.008*figSz(1) s4Pos2(2) 0.05*figSz(1) s4Pos2(4)]) '],'...
    '''FontSize'',' num2str(0.01*figSz(1)*txtScl) ','...
    '''FontWeight'',''Bold'','...
    '''HorizontalAlignment'',''Center'','...
    '''String'','''','...
    '''Callback'',@maxDispTextCallback_2,'...
    '''BackgroundColor'',[' num2str(get(findobj('Tag','3DviewerFig'),'color')) ']'];
s4t2 = initiateObject('s4Text2_2','uicontrol',str);
%S51 SLIDER
delete(findobj('Tag','s5Slider_1'));
s51 = uicontrol('Style','Slider','Position',s5Pos1,'Enable','off','Tag','s5Slider_1');
%S51 TEXT 1
str = ['''Style'',''PushButton'','...
    '''Units'',''Pixel'','...
    '''Position'',[' num2str([s4Pos1(1)-0.034*figSz(1) s5Pos1(2) 0.03*figSz(1) s5Pos1(4)]) '],'...
    '''FontSize'',' num2str(0.01*figSz(1)*txtScl) ','...
    '''HorizontalAlignment'',''Center'','...
    '''String'',''Min'','...
    '''Callback'',@minDispButtonCallback_1,'...
    '''BackgroundColor'',[' num2str(get(findobj('Tag','3DviewerFig'),'color')) ']'];
initiateObject('s5Button1_1','uicontrol',str);
%S51 TEXT 2
str = ['''Style'',''Edit'','...
    '''Units'',''Pixel'','...
    '''Position'',[' num2str([s5Pos1(1)+s5Pos1(3)+0.008*figSz(1) s5Pos1(2) 0.05*figSz(1) s5Pos1(4)]) '],'...
    '''FontSize'',' num2str(0.01*figSz(1)*txtScl) ','...
    '''FontWeight'',''Bold'','...
    '''HorizontalAlignment'',''Center'','...
    '''String'','''','...
    '''Callback'',@minDispTextCallback_1,'...
    '''BackgroundColor'',[' num2str(get(findobj('Tag','3DviewerFig'),'color')) ']'];
s5t1 = initiateObject('s5Text2_1','uicontrol',str);
%S52 SLIDER
delete(findobj('Tag','s5Slider_2'));
s52 = uicontrol('Style','Slider','Position',s5Pos2,'Enable','off','Tag','s5Slider_2');
%S52 TEXT 1
str = ['''Style'',''PushButton'','...
    '''Units'',''Pixel'','...
    '''Position'',[' num2str([s4Pos2(1)-0.034*figSz(1) s5Pos2(2) 0.03*figSz(1) s5Pos2(4)]) '],'...
    '''FontSize'',' num2str(0.01*figSz(1)*txtScl) ','...
    '''HorizontalAlignment'',''Center'','...
    '''String'',''Min'','...
    '''Callback'',@minDispButtonCallback_2 ,'...
    '''BackgroundColor'',[' num2str(get(findobj('Tag','3DviewerFig'),'color')) ']'];
initiateObject('s5Button1_2','uicontrol',str);
%S52 TEXT 2
str = ['''Style'',''Edit'','...
    '''Units'',''Pixel'','...
    '''Position'',[' num2str([s5Pos2(1)+s5Pos2(3)+0.008*figSz(1) s5Pos2(2) 0.05*figSz(1) s5Pos2(4)]) '],'...
    '''FontSize'',' num2str(0.01*figSz(1)*txtScl) ','...
    '''FontWeight'',''Bold'','...
    '''HorizontalAlignment'',''Center'','...
    '''String'','''','...
    '''Callback'',@minDispTextCallback_2,'...
    '''BackgroundColor'',[' num2str(get(findobj('Tag','3DviewerFig'),'color')) ']'];
s5t2 = initiateObject('s5Text2_2','uicontrol',str);

%add labels for each image
txtWdth = 0.12*figSz(1);
txtHght = 0.017*figSz(2);
pos1 = [s4Pos1(1)+(s4Pos1(3)-txtWdth)/2 s4Pos1(2)+s4Pos1(4)+(0.007*figSz(2)) txtWdth txtHght];
pos2 = [s4Pos2(1)+(s4Pos2(3)-txtWdth)/2 s4Pos2(2)+s4Pos2(4)+(0.007*figSz(2)) txtWdth txtHght];
im1str = '';
im2str = '';
str = ['''Style'',''Text'','...
    '''Units'',''Pixel'','...
    '''Position'',[' num2str(pos1) '],'...
    '''FontSize'',' num2str(0.006*figSz(1)*txtScl) ','...
    '''FontWeight'',''Bold'','...
    '''HorizontalAlignment'',''Center'','...
    '''String'',''' im1str ''','...
    '''BackgroundColor'',[' num2str(get(findobj('Tag','3DviewerFig'),'color')) ']'];
initiateObject('imLabel_1','uicontrol',str);
str = ['''Style'',''Text'','...
    '''Units'',''Pixel'','...
    '''Position'',[' num2str(pos2) '],'...
    '''FontSize'',' num2str(0.006*figSz(1)*txtScl) ','...
    '''FontWeight'',''Bold'','...
    '''HorizontalAlignment'',''Center'','...
    '''String'',''' im2str ''','...
    '''BackgroundColor'',[' num2str(get(findobj('Tag','3DviewerFig'),'color')) ']'];
initiateObject('imLabel_2','uicontrol',str);

end

function [s6,s6t] = formatSmoothingPanel(smpPos,s6Pos,txtScl)

if txtScl~=1 %high DPI scaling
    txtScl = txtScl * 0.7;
end

figSz = get(findobj('Tag','3DviewerFig'),'Position');
figSz = [figSz(3) figSz(4)];

%SLIDER PANEL
str = ['''Title'',''Gaussian Smoothing'','...
    '''Units'',''Pixel'','...
    '''Position'',[' num2str(smpPos) '],'...
    '''FontSize'',' num2str(0.01*figSz(1)*txtScl) ','...
    '''TitlePosition'',''CenterTop'','...
    '''BackgroundColor'',[' num2str(get(findobj('Tag','3DviewerFig'),'color')) ']'];
initiateObject('SmoothingPanel','uipanel',str);

%S6 SLIDER
delete(findobj('Tag','s6Slider'));
s6 = uicontrol('Style','Slider','Position',s6Pos,'Enable','off','Tag','s6Slider');
%S6 TEXT 1
str = ['''Style'',''Text'','...
    '''Units'',''Pixel'','...
    '''Position'',[' num2str([smpPos(1)+0.006*figSz(1) s6Pos(2) 0.047*figSz(1) s6Pos(4)]) '],'...
    '''FontSize'',' num2str(0.01*figSz(1)*txtScl) ','...
    '''HorizontalAlignment'',''Center'','...
    '''String'',''FWHM'','...
    '''BackgroundColor'',[' num2str(get(findobj('Tag','3DviewerFig'),'color')) ']'];
initiateObject('s6Text1','uicontrol',str);
%S6 TEXT 2
str = ['''Style'',''Edit'','...
    '''Units'',''Pixel'','...
    '''Position'',[' num2str([s6Pos(1)+s6Pos(3)+0.0065*figSz(1) s6Pos(2) 0.017*figSz(1) s6Pos(4)]) '],'...
    '''FontSize'',' num2str(0.01*figSz(1)*txtScl) ','...
    '''FontWeight'',''Bold'','...
    '''HorizontalAlignment'',''Center'','...
    '''String'','''','...
    '''Callback'',@smTextCallback,'...
    '''BackgroundColor'',[' num2str(get(findobj('Tag','3DviewerFig'),'color')) ']'];
s6t = initiateObject('s6Text2','uicontrol',str);
%S6 TEXT 2 UNITS (mm)
str = ['''Style'',''Text'','...
    '''Units'',''Pixel'','...
    '''Position'',[' num2str([s6Pos(1)+s6Pos(3)+0.025*figSz(1) s6Pos(2) 0.024*figSz(1) s6Pos(4)]) '],'...
    '''FontSize'',' num2str(0.01*figSz(1)*txtScl) ','...
    '''HorizontalAlignment'',''Center'','...
    '''String'',''mm'','...
    '''BackgroundColor'',[' num2str(get(findobj('Tag','3DviewerFig'),'color')) ']'];
initiateObject('s6Text2units','uicontrol',str);

end

function [s6] = formatFusedWeightPanel(fwpPos,s6Pos,txtScl)

if txtScl~=1 %high DPI scaling
    txtScl = txtScl * 0.7;
end

figSz = get(findobj('Tag','3DviewerFig'),'Position');
figSz = [figSz(3) figSz(4)];

%SLIDER PANEL
str = ['''Title'',''Image Fusion Weight'','...
    '''Units'',''Pixel'','...
    '''Position'',[' num2str(fwpPos) '],'...
    '''FontSize'',' num2str(0.01*figSz(1)*txtScl) ','...
    '''TitlePosition'',''CenterTop'','...
    '''BackgroundColor'',[' num2str(get(findobj('Tag','3DviewerFig'),'color')) ']'];
initiateObject('FusedWeightPanel','uipanel',str);

%S6 SLIDER
delete(findobj('Tag','s6Slider'))
%initiate with 50% alpha
s6 = uicontrol('Style','Slider','Position',s6Pos,...
    'Min',0,'Max',1,'Value',0.5,'Enable','off','Tag','s6Slider');
%S6 TEXT 1
str = ['''Style'',''Text'','...
    '''Units'',''Pixel'','...
    '''Position'',[' num2str([fwpPos(1)+0.006*figSz(1) s6Pos(2) 0.03*figSz(1) s6Pos(4)]) '],'...
    '''FontSize'',' num2str(0.01*figSz(1)*txtScl) ','...
    '''HorizontalAlignment'',''Center'','...
    '''String'','''','...
    '''BackgroundColor'',[' num2str(get(findobj('Tag','3DviewerFig'),'color')) ']'];
initiateObject('s6Text_1','uicontrol',str);
%S6 TEXT 2
str = ['''Style'',''Text'','...
    '''Units'',''Pixel'','...
    '''Position'',[' num2str([s6Pos(1)+s6Pos(3)+0.002*figSz(1) s6Pos(2) 0.03*figSz(1) s6Pos(4)]) '],'...
    '''FontSize'',' num2str(0.01*figSz(1)*txtScl) ','...
    '''HorizontalAlignment'',''Center'','...
    '''String'','''','...
    '''BackgroundColor'',[' num2str(get(findobj('Tag','3DviewerFig'),'color')) ']'];
initiateObject('s6Text_2','uicontrol',str);

end

function formatCMPanel(cmPos,txtScl,varargin)

%varargin is used for fused images to manually define the intance tags
if ~isempty(varargin)
    pTag = varargin{1};
    b1Tag = varargin{2};
    b2Tag = varargin{3};
    b3Tag = varargin{4};
else
    pTag = 'CMPanel';
    b1Tag = 'b1';
    b2Tag = 'b2';
    b3Tag = 'b3';
end

figSz = get(findobj('Tag','3DviewerFig'),'Position');
figSz = [figSz(3) figSz(4)];
bSz = [(-1.4e-5*figSz(1)+0.055)*figSz(1)/txtScl 0.03*figSz(2)];

if txtScl~=1 %high DPI scaling
    txtScl = txtScl * 0.7;
end

cmTxt = 'gray';
%ACCOUNT FOR PREVIOUSLY SET COLORMAP
test = findobj('Tag',pTag);
if ~isempty(test)
    test = get(test,'Children');
    if ~isempty(test)
        for n=1:numel(test)
            if get(test(n),'Value')==1
                cmTxt = lower(get(test(n),'String'));
            end
        end
    end
end
str = ['''Title'','''','...
    '''Units'',''Pixel'','...
    '''Position'',[' num2str(cmPos) '],'...
    '''FontSize'',' num2str(0.01*figSz(1)*txtScl) ','...
    '''TitlePosition'',''CenterTop'','...
    '''BorderType'',''None'','...
    '''SelectionChangeFcn'',@CMselect,'...
    '''BackgroundColor'',[' num2str(get(findobj('Tag','3DviewerFig'),'color')) ']'];
bg = initiateObject(pTag,'uibuttongroup',str);

str = ['''Style'',''Radiobutton'','...
    '''Units'',''Pixel'','...
    '''Position'',[' num2str([cmPos(3)/2-(bSz(1)/2) 0.07*figSz(2) bSz]) '],'...
    '''FontSize'',' num2str(0.006*figSz(1)*txtScl) ','...
    '''String'',''Gray'','...
    '''BackgroundColor'',[' num2str(get(findobj('Tag','3DviewerFig'),'color')) ']'];
initiateObject(b1Tag,'uicontrol',str,bg);
str = ['''Style'',''Radiobutton'','...
    '''Units'',''Pixel'','...
    '''Position'',[' num2str([cmPos(3)/2-(bSz(1)/2) 0.04*figSz(2) bSz]) '],'...
    '''FontSize'',' num2str(0.006*figSz(1)*txtScl) ','...
    '''String'',''Hot'','...
    '''BackgroundColor'',[' num2str(get(findobj('Tag','3DviewerFig'),'color')) ']'];
initiateObject(b2Tag,'uicontrol',str,bg);
str = ['''Style'',''Radiobutton'','...
    '''Units'',''Pixel'','...
    '''Position'',[' num2str([cmPos(3)/2-(bSz(1)/2) 0.01*figSz(2) bSz]) '],'...
    '''FontSize'',' num2str(0.006*figSz(1)*txtScl) ','...
    '''String'',''PET'','...
    '''BackgroundColor'',[' num2str(get(findobj('Tag','3DviewerFig'),'color')) ']'];
initiateObject(b3Tag,'uicontrol',str,bg);

test = get(bg,'Children');
set(test(strcmpi(get(test,'String'),cmTxt)),'Value',1)

if strcmpi(cmTxt,'gray'), setCM(1-gray); end

end

function formatCTPanel(ctPos,txtScl,varargin)

%varargin is used for fused images to manually define the intance tags
if ~isempty(varargin)
    pTag = varargin{1};
    b1Tag = varargin{2};
    b2Tag = varargin{3};
    b3Tag = varargin{4};
else
    pTag = 'CMPanel';
    b1Tag = 'b1';
    b2Tag = 'b2';
    b3Tag = 'b3';
end

figSz = get(findobj('Tag','3DviewerFig'),'Position');
figSz = [figSz(3) figSz(4)];
bSz = [(-1.4e-5*figSz(1)+0.055)*figSz(1)/txtScl 0.03*figSz(2)];

if txtScl~=1 %high DPI scaling
    txtScl = txtScl * 0.7;
end

cmTxt = 'body';
%ACCOUNT FOR PREVIOUSLY SET COLORMAP
test = findobj('Tag',pTag);
if ~isempty(test)
    test = get(test,'Children');
    if ~isempty(test)
        for n=1:numel(test)
            if get(test(n),'Value')==1
                cmTxt = lower(get(test(n),'String'));
            end
        end
    end
end

str = ['''Title'','''','...
    '''Units'',''Pixel'','...
    '''Position'',[' num2str(ctPos) '],'...
    '''FontSize'',' num2str(0.01*figSz(1)*txtScl) ','...
    '''TitlePosition'',''CenterTop'','...
    '''BorderType'',''None'','...
    '''SelectionChangeFcn'',@CTselect,'...
    '''BackgroundColor'',[' num2str(get(findobj('Tag','3DviewerFig'),'color')) ']'];
bg = initiateObject(pTag,'uibuttongroup',str);

str = ['''Style'',''Radiobutton'','...
    '''Units'',''Pixel'','...
    '''Position'',[' num2str([ctPos(3)/2-(bSz(1)/2) 0.07*figSz(2) bSz]) '],'...
    '''FontSize'',' num2str(0.006*figSz(1)*txtScl) ','...
    '''String'',''Body'','...
    '''BackgroundColor'',[' num2str(get(findobj('Tag','3DviewerFig'),'color')) ']'];
initiateObject(b1Tag,'uicontrol',str,bg);
str = ['''Style'',''Radiobutton'','...
    '''Units'',''Pixel'','...
    '''Position'',[' num2str([ctPos(3)/2-(bSz(1)/2) 0.04*figSz(2) bSz]) '],'...
    '''FontSize'',' num2str(0.006*figSz(1)*txtScl) ','...
    '''String'',''Lung'','...
    '''BackgroundColor'',[' num2str(get(findobj('Tag','3DviewerFig'),'color')) ']'];
initiateObject(b2Tag,'uicontrol',str,bg);
str = ['''Style'',''Radiobutton'','...
    '''Units'',''Pixel'','...
    '''Position'',[' num2str([ctPos(3)/2-(bSz(1)/2) 0.01*figSz(2) bSz]) '],'...
    '''FontSize'',' num2str(0.006*figSz(1)*txtScl) ','...
    '''String'',''Bone'','...
    '''BackgroundColor'',[' num2str(get(findobj('Tag','3DviewerFig'),'color')) ']'];
initiateObject(b3Tag,'uicontrol',str,bg);

test = get(bg,'Children');
set(test(strcmpi({get(test,'String')},cmTxt)),'Value',1)

end

function makeActive

global handles

%FIGURES AND AXES
if ~handles.display.images.fusedImages
    set(handles.figure.f1,'WindowButtonMotionFcn',@figMotionFcn,'ButtonDownFcn',@toggleFcn,...
        'WindowButtonUpFcn',@stopDragging,'WindowScrollWheelFcn',@zoomFcn,...
        'ResizeFcn',@figResize,'Interruptible','Off','BusyAction','Queue')
else
    set(handles.figure.f1,'WindowButtonMotionFcn',@figMotionFcnFused,'ButtonDownFcn',@toggleFcn,...
        'WindowButtonUpFcn',@stopDragging,'WindowScrollWheelFcn',@zoomFcn,...
        'ResizeFcn',@figResize,'Interruptible','Off','BusyAction','Queue')
end
set(handles.axes.a1,'ButtonDownFcn',@a1ButtonDownFcn,'Visible','on')
set(get(handles.axes.a1,'Title'),'String',handles.axes.a1title)
set(findobj('Tag','IOPt1'),'Visible','on')
set(findobj('Tag','IOPt2'),'Visible','on')
set(handles.axes.a2,'ButtonDownFcn',@a2ButtonDownFcn,'Visible','on')
set(get(handles.axes.a2,'Title'),'String',handles.axes.a2title)
set(findobj('Tag','IOPt3'),'Visible','on')
set(findobj('Tag','IOPt4'),'Visible','on')
set(handles.axes.a3,'ButtonDownFcn',@a3ButtonDownFcn,'Visible','on')
set(get(handles.axes.a3,'Title'),'String',handles.axes.a3title)
set(findobj('Tag','IOPt5'),'Visible','on')
set(findobj('Tag','IOPt6'),'Visible','on')
%AXES
% set(findobj('Tag','axes1'),'Visible','on')
% set(findobj('Tag','axes2'),'Visible','on')
% set(findobj('Tag','axes3'),'Visible','on')
%IMAGES
set(findobj('Tag','a1im'),'Visible','on')
set(findobj('Tag','a2im'),'Visible','on')
set(findobj('Tag','a3im'),'Visible','on')
set(findobj('Tag','a1im_1'),'Visible','on')
set(findobj('Tag','a1im_2'),'Visible','on')
set(findobj('Tag','a2im_1'),'Visible','on')
set(findobj('Tag','a2im_2'),'Visible','on')
set(findobj('Tag','a3im_1'),'Visible','on')
set(findobj('Tag','a3im_2'),'Visible','on')
%LOCATION ARROWS
set(findall(0,'Tag','a1arrowH1'),'Visible','on')
set(findall(0,'Tag','a1arrowH2'),'Visible','on')
set(findall(0,'Tag','a1arrowV1'),'Visible','on')
set(findall(0,'Tag','a1arrowV2'),'Visible','on')
set(findall(0,'Tag','a2arrowH1'),'Visible','on')
set(findall(0,'Tag','a2arrowH2'),'Visible','on')
set(findall(0,'Tag','a2arrowV1'),'Visible','on')
set(findall(0,'Tag','a2arrowV2'),'Visible','on')
set(findall(0,'Tag','a3arrowH1'),'Visible','on')
set(findall(0,'Tag','a3arrowH2'),'Visible','on')
set(findall(0,'Tag','a3arrowV1'),'Visible','on')
set(findall(0,'Tag','a3arrowV2'),'Visible','on')
%NAVIGATION SLIDERS
if handles.volume.imSz3>1
    set(handles.display.sliders.s1,'Enable','on','Callback',@releaseSlider)
    if exist('handles.display.sliders.s1listener','var')
        delete(handles.display.sliders.s1listener)
    end
    handles.display.sliders.s1listener = ...
        addlistener(handles.display.sliders.s1,'ContinuousValueChange',@a1SliderListener);
end
set(handles.display.sliders.s1t,'Enable','on','Visible','on')
if handles.volume.imSz2>1
    set(handles.display.sliders.s2,'Enable','on','Callback',@releaseSlider)
    if exist('handles.display.sliders.s2listener','var')
        delete(handles.display.sliders.s2listener)
    end
    handles.display.sliders.s2listener = ...
        addlistener(handles.display.sliders.s2,'ContinuousValueChange',@a2SliderListener);
end
set(handles.display.sliders.s2t,'Enable','on','Visible','on')
if handles.volume.imSz1>1
    set(handles.display.sliders.s3,'Enable','on','Callback',@releaseSlider)
    if exist('handles.display.sliders.s3listener','var')
        delete(handles.display.sliders.s3listener)
    end
    handles.display.sliders.s3listener = ...
        addlistener(handles.display.sliders.s3,'ContinuousValueChange',@a3SliderListener);
end
set(handles.display.sliders.s3t,'Enable','on','Visible','on')
if ~handles.display.images.fusedImages
    %DISPLAY SLIDERS
    set(handles.display.sliders.s4,'Enable','on'); %,'Callback',@maxDispSliderCallback)
    if exist('handles.display.sliders.s4listener','var')
        delete(handles.display.sliders.s4listener)
    end
    handles.display.sliders.s4listener = ...
        addlistener(handles.display.sliders.s4,'ContinuousValueChange',@maxDispSliderCallback);
    set(handles.display.sliders.s4t,'Enable','on','Visible','on')
    set(handles.display.sliders.s5,'Enable','on'); %,'Callback',@minDispSliderCallback)
    if exist('handles.display.sliders.s5listener','var')
        delete(handles.display.sliders.s5listener)
    end
    handles.display.sliders.s5listener = ...
        addlistener(handles.display.sliders.s5,'ContinuousValueChange',@minDispSliderCallback);
    set(handles.display.sliders.s5t,'Enable','on','Visible','on')
    %SMOOTHING SLIDER
    set(handles.display.sliders.s6,'Enable','on','Callback',@smSliderCallback)
    if exist('handles.display.sliders.s6listener','var')
        delete(handles.display.sliders.s6listener)
    end
    handles.display.sliders.s6listener = ...
        addlistener(handles.display.sliders.s6,'ContinuousValueChange',@smSliderText);
    set(handles.display.sliders.s6t,'Enable','on','Visible','on')
    set(findobj('Tag','s4Button1'),'Enable','on')
    set(findobj('Tag','s5Button1'),'Enable','on')
    set(findobj('Tag','s6Text1'),'Enable','on')
    set(findobj('Tag','s6Text2'),'Enable','on')
    set(findobj('Tag','s6Text2units'),'Enable','on')
    %COLORMAP PANEL
    set(findobj('Tag','b1'),'Enable','on')
    set(findobj('Tag','b2'),'Enable','on')
    set(findobj('Tag','b3'),'Enable','on')
else
    %DISPLAY SLIDERS
    set(handles.display.sliders.s4_1,'Enable','on');
    if exist('handles.display.sliders.s4listener_1','var')
        delete(handles.display.sliders.s4listener_1)
    end
    handles.display.sliders.s4listener_1 = ...
        addlistener(handles.display.sliders.s4_1,'ContinuousValueChange',@maxDispSliderCallback_1);
    set(handles.display.sliders.s4t_1,'Enable','on','Visible','on')
    set(findobj('Tag','s4Button1_1'),'Enable','on')
    set(handles.display.sliders.s5_1,'Enable','on');
    if exist('handles.display.sliders.s5listener_1','var')
        delete(handles.display.sliders.s5listener_1)
    end
    handles.display.sliders.s5listener_1 = ...
        addlistener(handles.display.sliders.s5_1,'ContinuousValueChange',@minDispSliderCallback_1);
    set(handles.display.sliders.s5t_1,'Enable','on','Visible','on')
    set(findobj('Tag','s5Button1_1'),'Enable','on')
    
    set(handles.display.sliders.s4_2,'Enable','on');
    if exist('handles.display.sliders.s4listener_2','var')
        delete(handles.display.sliders.s4listener_2)
    end
    handles.display.sliders.s4listener_2 = ...
        addlistener(handles.display.sliders.s4_2,'ContinuousValueChange',@maxDispSliderCallback_2);
    set(handles.display.sliders.s4t_2,'Enable','on','Visible','on')
    set(findobj('Tag','s4Button1_2'),'Enable','on')
    set(handles.display.sliders.s5_2,'Enable','on');
    if exist('handles.display.sliders.s5listener_2','var')
        delete(handles.display.sliders.s5listener_2)
    end
    handles.display.sliders.s5listener_2 = ...
        addlistener(handles.display.sliders.s5_2,'ContinuousValueChange',@minDispSliderCallback_2);
    set(handles.display.sliders.s5t_2,'Enable','on','Visible','on')
    set(findobj('Tag','s5Button1_2'),'Enable','on')
    
    set(findobj('Tag','imLabel_1'),'Enable','on')
    set(findobj('Tag','imLabel_2'),'Enable','on')
    
    %FUSION WEIGHT SLIDER
    set(handles.display.sliders.s6_1,'Enable','on','Callback',@fwSliderCallback)
    if exist('handles.display.sliders.s6listener_1','var')
        delete(handles.display.sliders.s6listener_1)
    end
    handles.display.sliders.s6listener_1 = ...
        addlistener(handles.display.sliders.s6_1,'ContinuousValueChange',@fwSliderCallback);
    set(findobj('Tag','s6Text_1'),'Enable','on')
    set(findobj('Tag','s6Text_2'),'Enable','on')
    %VOXEL TEXT
    set(findobj('Tag','imVoxLabel_1'),'Enable','on')
    set(findobj('Tag','imVoxLabel_2'),'Enable','on')
end

%CM COLORMAP OR CT PRESET PANEL
set(findobj('Tag','b1_1'),'Enable','on')
set(findobj('Tag','b2_1'),'Enable','on')
set(findobj('Tag','b3_1'),'Enable','on')
set(findobj('Tag','b1_2'),'Enable','on')
set(findobj('Tag','b2_2'),'Enable','on')
set(findobj('Tag','b3_2'),'Enable','on')

set(findobj('Tag','LocText1x'),'Visible','on')
set(findobj('Tag','LocText1y'),'Visible','on')
set(findobj('Tag','LocText1z'),'Visible','on')
set(findobj('Tag','LocText2x'),'Visible','on')
set(findobj('Tag','LocText2y'),'Visible','on')
set(findobj('Tag','LocText2z'),'Visible','on')

%VOLUME STAMP TEXT
set(findobj('Tag','StampText1'),'Visible','on')
set(findobj('Tag','StampText2'),'Visible','on')

if isfield(handles,'drawBtn')
    handles.mask.draw = false;
    set(handles.drawBtn,'Enable','on','Callback',@drawButtonFcn,'BusyAction','Cancel')
end
if isfield(handles,'doneBtn')
    set(handles.doneBtn,'Enable','on','Callback','uiresume','BusyAction','Cancel')
end

end

function makeInactive

%FIGURES AND AXES
set(findobj('Tag','3DviewerFig'),'WindowButtonMotionFcn','',...
    'ResizeFcn','','ButtonDownFcn','','WindowButtonUpFcn','')
test = findobj('Tag','axes1');
if ~isempty(test)
    set(test,'ButtonDownFcn','','XColor',[0 0 0],'YColor',[0 0 0]);
    set(get(test,'Title'),'String','')
end
set(findobj('Tag','IOPt1'),'Visible','off')
set(findobj('Tag','IOPt2'),'Visible','off')
test = findobj('Tag','axes2');
if ~isempty(test)
    set(test,'ButtonDownFcn','','XColor',[0 0 0],'YColor',[0 0 0]);
    set(get(test,'Title'),'String','')
end
set(findobj('Tag','IOPt3'),'Visible','off')
set(findobj('Tag','IOPt4'),'Visible','off')
test = findobj('Tag','axes3');
if ~isempty(test)
    set(test,'ButtonDownFcn','','XColor',[0 0 0],'YColor',[0 0 0]);
    set(get(test,'Title'),'String','')
end
set(findobj('Tag','IOPt5'),'Visible','off')
set(findobj('Tag','IOPt6'),'Visible','off')
%AXES
% set(findobj('Tag','axes1'),'Visible','off')
% set(findobj('Tag','axes2'),'Visible','off')
% set(findobj('Tag','axes3'),'Visible','off')
%IMAGES
set(findobj('Tag','a1im'),'Visible','off')
set(findobj('Tag','a2im'),'Visible','off')
set(findobj('Tag','a3im'),'Visible','off')
set(findobj('Tag','a1im_1'),'Visible','off')
set(findobj('Tag','a1im_2'),'Visible','off')
set(findobj('Tag','a2im_1'),'Visible','off')
set(findobj('Tag','a2im_2'),'Visible','off')
set(findobj('Tag','a3im_1'),'Visible','off')
set(findobj('Tag','a3im_2'),'Visible','off')
%LOCATION ARROWS
set(findall(0,'Tag','a1arrowH1'),'Visible','off')
set(findall(0,'Tag','a1arrowH2'),'Visible','off')
set(findall(0,'Tag','a1arrowV1'),'Visible','off')
set(findall(0,'Tag','a1arrowV2'),'Visible','off')
set(findall(0,'Tag','a2arrowH1'),'Visible','off')
set(findall(0,'Tag','a2arrowH2'),'Visible','off')
set(findall(0,'Tag','a2arrowV1'),'Visible','off')
set(findall(0,'Tag','a2arrowV2'),'Visible','off')
set(findall(0,'Tag','a3arrowH1'),'Visible','off')
set(findall(0,'Tag','a3arrowH2'),'Visible','off')
set(findall(0,'Tag','a3arrowV1'),'Visible','off')
set(findall(0,'Tag','a3arrowV2'),'Visible','off')
%NAVIGATION SLIDERS
set(findobj('Tag','s1Slider'),'Callback','','Enable','off')
set(findobj('Tag','s1Text'),'String','','Enable','off')
set(findobj('Tag','s2Slider'),'Callback','','Enable','off')
set(findobj('Tag','s2Text'),'String','','Enable','off')
set(findobj('Tag','s3Slider'),'Callback','','Enable','off')
set(findobj('Tag','s3Text'),'String','','Enable','off')
%DISPLAY SLIDERS
set(findobj('Tag','s4Slider'),'Callback','','Enable','off')
set(findobj('Tag','s4Button1'),'Enable','off')
set(findobj('Tag','s4Text2'),'String','','Enable','off')
set(findobj('Tag','s5Slider'),'Callback','','Enable','off')
set(findobj('Tag','s5Button1'),'Enable','off')
set(findobj('Tag','s5Text2'),'String','','Enable','off')
set(findobj('Tag','s4Slider_1'),'Callback','','Enable','off')
set(findobj('Tag','s4Button1_1'),'Enable','off')
set(findobj('Tag','s4Text2_1'),'String','','Enable','off')
set(findobj('Tag','s5Slider_1'),'Callback','','Enable','off')
set(findobj('Tag','s5Button1_1'),'Enable','off')
set(findobj('Tag','s5Text2_1'),'String','','Enable','off')
set(findobj('Tag','s4Slider_2'),'Callback','','Enable','off')
set(findobj('Tag','s4Button1_2'),'Enable','off')
set(findobj('Tag','s4Text2_2'),'String','','Enable','off')
set(findobj('Tag','s5Slider_2'),'Callback','','Enable','off')
set(findobj('Tag','s5Button1_2'),'Enable','off')
set(findobj('Tag','s5Text2_2'),'String','','Enable','off')
set(findobj('Tag','imLabel_1'),'String','','Enable','off')
set(findobj('Tag','imLabel_2'),'String','','Enable','off')
%SMOOTHING SLIDER
set(findobj('Tag','s6Slider'),'Callback','','Enable','off')
set(findobj('Tag','s6Text1'),'Enable','off')
set(findobj('Tag','s6Text2'),'String','','Enable','off')
set(findobj('Tag','s6Text2units'),'Enable','off')
set(findobj('Tag','s6Text_1'),'Enable','off')
set(findobj('Tag','s6Text_2'),'Enable','off')
%COLORMAP PANEL
set(findobj('Tag','b1'),'Enable','off')
set(findobj('Tag','b2'),'Enable','off')
set(findobj('Tag','b3'),'Enable','off')
set(findobj('Tag','b1_1'),'Enable','off')
set(findobj('Tag','b2_1'),'Enable','off')
set(findobj('Tag','b3_1'),'Enable','off')
set(findobj('Tag','b1_2'),'Enable','off')
set(findobj('Tag','b2_2'),'Enable','off')
set(findobj('Tag','b3_2'),'Enable','off')
%VOXEL AND COORDINATE TEXT
set(findobj('Tag','VoxText'),'Visible','off')
set(findobj('Tag','VoxText_1'),'Visible','off')
set(findobj('Tag','VoxText_2'),'Visible','off')
set(findobj('Tag','imVoxLabel_1'),'Enable','off')
set(findobj('Tag','imVoxLabel_2'),'Enable','off')
set(findobj('Tag','CoordinateText'),'Visible','off')
set(findobj('Tag','LocText1x'),'Visible','off')
set(findobj('Tag','LocText1y'),'Visible','off')
set(findobj('Tag','LocText1z'),'Visible','off')
set(findobj('Tag','LocText2x'),'Visible','off')
set(findobj('Tag','LocText2y'),'Visible','off')
set(findobj('Tag','LocText2z'),'Visible','off')
%VOLUME STAMP TEXT
set(findobj('Tag','StampText1'),'Visible','off')
set(findobj('Tag','StampText2'),'Visible','off')

set(findobj('Tag','maskPatch1'),'Visible','off')
set(findobj('Tag','maskPatch2'),'Visible','off')
set(findobj('Tag','maskPatch3'),'Visible','off')

end

function makeFigScalable

global handles

figPos = get(handles.figure.f1,'Position');

temp = findobj(handles.figure.f1,'-not','type','Root','-not','type','Figure',...
    '-not','type','Line','-not','type','Patch','-not','type','Image');
set(temp,'Units','Normalized')

%CONVERT SLIDER TEXTS TO NORMALIZED FIGURE COORDINATES
handles.display.sliders.s1tLocs = handles.display.sliders.s1tLocs/figPos(3);
handles.display.sliders.s2tLocs = handles.display.sliders.s2tLocs/figPos(3);
handles.display.sliders.s3tLocs = handles.display.sliders.s3tLocs/figPos(3);
handles.display.sliders.st_y = handles.display.sliders.st_y/figPos(4);
handles.display.sliders.stSz = handles.display.sliders.stSz./[figPos(3) figPos(4)];

set(handles.display.sliders.s1t,'Units','Normalized');
set(handles.display.sliders.s2t,'Units','Normalized');
set(handles.display.sliders.s3t,'Units','Normalized');

handles.axes.a1Pos = get(handles.axes.a1,'Position')+[-0.001 -0.001 0.002 0.002];
handles.axes.a2Pos = get(handles.axes.a2,'Position')+[-0.001 -0.001 0.002 0.002];
handles.axes.a3Pos = get(handles.axes.a3,'Position')+[-0.001 -0.001 0.002 0.002];
handles.display.sliders.s1Pos = get(handles.display.sliders.s1,'Position');
handles.display.sliders.s2Pos = get(handles.display.sliders.s2,'Position');
handles.display.sliders.s3Pos = get(handles.display.sliders.s3,'Position');
if ~handles.display.images.fusedImages
    handles.display.sliders.s4Pos = get(handles.display.sliders.s4,'Position');
    handles.display.sliders.s5Pos = get(handles.display.sliders.s5,'Position');
    handles.display.sliders.s6Pos = get(handles.display.sliders.s6,'Position');
else
    handles.display.sliders.s4Pos_1 = get(handles.display.sliders.s4_1,'Position');
    handles.display.sliders.s4Pos_2 = get(handles.display.sliders.s4_2,'Position');
    handles.display.sliders.s5Pos_1 = get(handles.display.sliders.s5_1,'Position');
    handles.display.sliders.s5Pos_2 = get(handles.display.sliders.s5_2,'Position');
    handles.display.sliders.s6Pos_1 = get(handles.display.sliders.s6_1,'Position');
end

end

function figResize(hObj,~)

%this function is unused because the resize commands are commented out

set(hObj,'Interruptible','On')

try %this avoids errors reported on some platforms
    
    global handles
    
    oldPos = handles.figure.figPos;
    newPos = get(handles.figure.f1,'Position');
    
    sizeChg = abs(newPos(3:4)-oldPos(3:4)); % Change in figure size: [widthchg heightchg]
    if sizeChg(1) >= sizeChg(2)
        % Width change is larger, so change height to keep figure aspect ratio constant
        newPos(3:4) = [newPos(3) newPos(3)/handles.figure.aspectRatio];
    else
        % Height change is larger, so change width to keep figure aspect ratio constant
        newPos(3:4) = [newPos(4)*handles.figure.aspectRatio newPos(4)];
    end
    
    figSz = newPos(3:4);
    
    fs1 = 0.01 * figSz(1) * handles.figure.DPIscl;
    if handles.figure.DPIscl~=1
        fs1 = fs1 * 0.7;
    end
    fs5 = 0.005*figSz(1)*handles.figure.DPIscl;
    fs6 = 0.006 * figSz(1) * handles.figure.DPIscl;
    if handles.figure.DPIscl~=1
        fs6 = fs6 * 0.7;
    end
    fs7 = 0.007*figSz(1)*handles.figure.DPIscl;
    fs117 = 0.0117*figSz(1)*handles.figure.DPIscl;
    if handles.figure.DPIscl~=1
        fs117 = fs117 * 0.7;
    end
    
    title(handles.axes.a1,get(get(handles.axes.a1,'Title'),'String'),'FontSize',fs1);
    title(handles.axes.a2,get(get(handles.axes.a2,'Title'),'String'),'FontSize',fs1);
    title(handles.axes.a3,get(get(handles.axes.a3,'Title'),'String'),'FontSize',fs1);
    set(findobj('Tag','InstructionsText'),'FontSize',fs7)
    set(findobj('Tag','VoxPanel'),'FontSize',fs1)
    set(findobj('Tag','VoxText'),'FontSize',fs1)
    set(findobj('Tag','VoxText_1'),'FontSize',fs1)
    set(findobj('Tag','VoxText_2'),'FontSize',fs1)
    set(findobj('Tag','imVoxLabel_1'),'FontSize',fs1)
    set(findobj('Tag','imVoxLabel_2'),'FontSize',fs1)
    set(findobj('Tag','CoordinatePanel'),'FontSize',fs1)
    set(findobj('Tag','CoordinateText'),'FontSize',fs1)
    set(findobj('Tag','LocText1x'),'FontSize',fs7)
    set(findobj('Tag','LocText1y'),'FontSize',fs7)
    set(findobj('Tag','LocText1z'),'FontSize',fs7)
    set(findobj('Tag','LocText2x'),'FontSize',fs7)
    set(findobj('Tag','LocText2y'),'FontSize',fs7)
    set(findobj('Tag','LocText2z'),'FontSize',fs7)
    set(findobj('Tag','DisplayPanel'),'FontSize',fs1)
    set(findobj('Tag','s1Text'),'FontSize',fs7)
    set(findobj('Tag','s2Text'),'FontSize',fs7)
    set(findobj('Tag','s3Text'),'FontSize',fs7)
    set(findobj('Tag','s4Button1'),'FontSize',fs1)
    set(findobj('Tag','s4Text2'),'FontSize',fs1)
    set(findobj('Tag','s5Button1'),'FontSize',fs1)
    set(findobj('Tag','s5Text2'),'FontSize',fs1)
    set(findobj('Tag','s4Button1_1'),'FontSize',fs1)
    set(findobj('Tag','s4Text2_1'),'FontSize',fs1)
    set(findobj('Tag','s5Button1_1'),'FontSize',fs1)
    set(findobj('Tag','s5Text2_1'),'FontSize',fs1)
    set(findobj('Tag','s4Button1_2'),'FontSize',fs1)
    set(findobj('Tag','s4Text2_2'),'FontSize',fs1)
    set(findobj('Tag','s5Button1_2'),'FontSize',fs1)
    set(findobj('Tag','s5Text2_2'),'FontSize',fs1)
    set(findobj('Tag','imLabel_1'),'FontSize',fs6)
    set(findobj('Tag','imLabel_2'),'FontSize',fs6)
    set(findobj('Tag','SmoothingPanel'),'FontSize',fs1)
    set(findobj('Tag','FusedWeightPanel'),'FontSize',fs1)
    set(findobj('Tag','s6Text1'),'FontSize',fs1)
    set(findobj('Tag','s6Text2'),'FontSize',fs1)
    set(findobj('Tag','s6Text_1'),'FontSize',fs1)
    set(findobj('Tag','s6Text_2'),'FontSize',fs1)
    set(findobj('Tag','CMPanel'),'FontSize',fs1)
    set(findobj('Tag','b1'),'FontSize',fs6)
    set(findobj('Tag','b2'),'FontSize',fs6)
    set(findobj('Tag','b3'),'FontSize',fs6)
    set(findobj('Tag','b1_1'),'FontSize',fs6)
    set(findobj('Tag','b2_1'),'FontSize',fs6)
    set(findobj('Tag','b3_1'),'FontSize',fs6)
    set(findobj('Tag','b1_2'),'FontSize',fs6)
    set(findobj('Tag','b2_2'),'FontSize',fs6)
    set(findobj('Tag','b3_2'),'FontSize',fs6)
    set(findobj('Tag','IOPt1'),'FontSize',fs6)
    set(findobj('Tag','IOPt2'),'FontSize',fs6)
    set(findobj('Tag','IOPt3'),'FontSize',fs6)
    set(findobj('Tag','IOPt4'),'FontSize',fs6)
    set(findobj('Tag','IOPt5'),'FontSize',fs6)
    set(findobj('Tag','IOPt6'),'FontSize',fs6)
    set(findobj('Tag','doneBtn'),'FontSize',fs1)
    set(findobj('Tag','StampText1'),'FontSize',fs5)
    set(findobj('Tag','StampText2'),'FontSize',fs5)
    
    %NOW ARROWS
    headLength = 0.0062*figSz(1)*handles.figure.DPIscl;
    headWidth = 0.0047*figSz(1)*handles.figure.DPIscl;
    set(findall(0,'Tag','a1arrowH1'),'HeadLength',headLength,'HeadWidth',headWidth)
    set(findall(0,'Tag','a1arrowH2'),'HeadLength',headLength,'HeadWidth',headWidth)
    set(findall(0,'Tag','a1arrowV1'),'HeadLength',headLength,'HeadWidth',headWidth)
    set(findall(0,'Tag','a1arrowV2'),'HeadLength',headLength,'HeadWidth',headWidth)
    set(findall(0,'Tag','a2arrowH1'),'HeadLength',headLength,'HeadWidth',headWidth)
    set(findall(0,'Tag','a2arrowH2'),'HeadLength',headLength,'HeadWidth',headWidth)
    set(findall(0,'Tag','a2arrowV1'),'HeadLength',headLength,'HeadWidth',headWidth)
    set(findall(0,'Tag','a2arrowV2'),'HeadLength',headLength,'HeadWidth',headWidth)
    set(findall(0,'Tag','a3arrowH1'),'HeadLength',headLength,'HeadWidth',headWidth)
    set(findall(0,'Tag','a3arrowH2'),'HeadLength',headLength,'HeadWidth',headWidth)
    set(findall(0,'Tag','a3arrowV1'),'HeadLength',headLength,'HeadWidth',headWidth)
    set(findall(0,'Tag','a3arrowV2'),'HeadLength',headLength,'HeadWidth',headWidth)
    
    % set(hObj,'Position',newPos)
    % handles.figure.figPos = newPos;

catch
    
end

set(hObj,'Interruptible','Off')
% findfigs

end

function loadVolume(in)

global handles

figSz = get(handles.figure.f1,'Position');
figSz = [figSz(3) figSz(4)];

if ~isfield(in,'FOV1')
    in.FOV1 = in.voxDim1*size(in.vol,1);
end
if ~isfield(in,'FOV2')
    in.FOV2 = in.voxDim2*size(in.vol,2);
end
if ~isfield(in,'FOV3')
    in.FOV3 = in.voxDim3*size(in.vol,3);
end

if in.FOV1>in.FOV3
    scl = (handles.figure.ax3x-handles.figure.ax1x-2*handles.figure.axGap)/(2*in.FOV1);
elseif in.FOV2>in.FOV3
    scl = (handles.figure.ax3x-handles.figure.ax1x-2*handles.figure.axGap)/(2*in.FOV2);
else
    figPos = get(handles.figure.f1,'Position');
    scl = 0.5*figPos(4)/in.FOV3;
end
ax1Pos = [handles.figure.ax1x-scl*in.FOV1/2 handles.figure.ax_y-scl*in.FOV2/2 scl*in.FOV1 scl*in.FOV2];
ax2Pos = [handles.figure.ax2x-scl*in.FOV1/2 handles.figure.ax_y-scl*in.FOV3/2 scl*in.FOV1 scl*in.FOV3];
ax3Pos = [handles.figure.ax3x-scl*in.FOV2/2 handles.figure.ax_y-scl*in.FOV3/2 scl*in.FOV2 scl*in.FOV3];
while ax1Pos(1)+ax1Pos(3)>=ax2Pos(1)-handles.figure.axGap/2 || ...
        ax2Pos(1)+ax2Pos(3)>=ax3Pos(1)-handles.figure.axGap/2
    scl = scl-0.1;
    ax1Pos = [handles.figure.ax1x-scl*in.FOV1/2 handles.figure.ax_y-scl*in.FOV2/2 scl*in.FOV1 scl*in.FOV2];
    ax2Pos = [handles.figure.ax2x-scl*in.FOV1/2 handles.figure.ax_y-scl*in.FOV3/2 scl*in.FOV1 scl*in.FOV3];
    ax3Pos = [handles.figure.ax3x-scl*in.FOV2/2 handles.figure.ax_y-scl*in.FOV3/2 scl*in.FOV2 scl*in.FOV3];
end
[handles.axes.a1,handles.axes.a2,handles.axes.a3] = formatAxes(ax1Pos,ax2Pos,ax3Pos);

%use volume of absolute voxel locations
if isfield(in,'locX')
    handles.display.text.showLocs = true;
    handles.volume.locX = in.locX;
    handles.volume.locY = in.locY;
    handles.volume.locZ = in.locZ;
    %initiate text fields to display voxel locations
    [handles.display.text.locTxtX,...
        handles.display.text.locTxtY,...
        handles.display.text.locTxtZ] = formatLocPanel();
end

if ~handles.display.images.fusedImages
    handles.volume.vol0mm = in.vol;
    handles.volume.imSz1 = size(in.vol,1);
    handles.volume.imSz2 = size(in.vol,2);
    handles.volume.imSz3 = size(in.vol,3);
    handles.volume.voxDimX = in.FOV1/handles.volume.imSz1;
    handles.volume.voxDimY = in.FOV2/handles.volume.imSz2;
    handles.volume.voxDimZ = in.FOV3/handles.volume.imSz3;
    %DEFINE X,Y,Z, START WITH CENTER SLICES
    handles.display.images.x = round(handles.volume.imSz1/2);
    handles.display.images.y = round(handles.volume.imSz2/2);
    handles.display.images.z = round(handles.volume.imSz3/2);
    %DEFINE SMOOTHING SLIDERS
    maxVal = 6; %mm FWHM
    if get(handles.display.sliders.s6,'Value')~=0
        %     val = round(get(handles.display.sliders.s6,'Value'));
        val = 0;
    else
        val = 0;
    end
    set(handles.display.sliders.s6,'Min',0,'Max',maxVal,'Value',val,'SliderStep',[1/maxVal 2/maxVal])
    set(handles.display.sliders.s6t,'String',num2str(val))
    %ACCOUNT FOR ASYMMETRIC VOXELS
    siz = val./[handles.volume.voxDimX handles.volume.voxDimY handles.volume.voxDimZ];
    %SMOOTH VOLUME
    handles.volume.vol = GaussianSmooth(handles.volume.vol0mm, siz);
    %LOAD VOLUME INTO WINDOWS
    handles.display.images.h1 = setVolume(findobj('Tag','a1im'),handles.axes.a1);
    handles.display.images.h2 = setVolume(findobj('Tag','a2im'),handles.axes.a2);
    handles.display.images.h3 = setVolume(findobj('Tag','a3im'),handles.axes.a3);
    %clim
    cLimMin = min(handles.volume.vol(:));
    cLimMax = max(handles.volume.vol(:));
else
    handles.volume.vol_1 = in.vol_1;
    handles.volume.vol_2 = in.vol_2;
    %volume matrix dimensions are the same
    handles.volume.imSz1 = size(in.vol_1,1);
    handles.volume.imSz2 = size(in.vol_1,2);
    handles.volume.imSz3 = size(in.vol_1,3);
    handles.volume.voxDimX = in.FOV1/handles.volume.imSz1;
    handles.volume.voxDimY = in.FOV2/handles.volume.imSz2;
    handles.volume.voxDimZ = in.FOV3/handles.volume.imSz3;
    %DEFINE X,Y,Z, START WITH CENTER SLICES
    handles.display.images.x = round(handles.volume.imSz1/2);
    handles.display.images.y = round(handles.volume.imSz2/2);
    handles.display.images.z = round(handles.volume.imSz3/2);
    %clim
    cLimMin_1 = min(handles.volume.vol_1(:));
    cLimMax_1 = max(handles.volume.vol_1(:));
    cLimMin_2 = min(handles.volume.vol_2(:));
    cLimMax_2 = max(handles.volume.vol_2(:));
    %load volumes
    [handles.display.images.h1_1,handles.display.images.h1_2] = ...
        setFusedVolume(findobj('Tag','a1im_1'),findobj('Tag','a1im_2'),handles.axes.a1,...
        [cLimMin_1 cLimMax_1],[cLimMin_2 cLimMax_2],...
        handles.display.images.cm_1,handles.display.images.cm_2);
    [handles.display.images.h2_1,handles.display.images.h2_2] = ...
        setFusedVolume(findobj('Tag','a2im_1'),findobj('Tag','a2im_2'),handles.axes.a2,...
        [cLimMin_1 cLimMax_1],[cLimMin_2 cLimMax_2],...
        handles.display.images.cm_1,handles.display.images.cm_2);
    [handles.display.images.h3_1,handles.display.images.h3_2] = ...
        setFusedVolume(findobj('Tag','a3im_1'),findobj('Tag','a3im_2'),handles.axes.a3,...
        [cLimMin_1 cLimMax_1],[cLimMin_2 cLimMax_2],...
        handles.display.images.cm_1,handles.display.images.cm_2);
end

handles.display.images.xLimits = [1 handles.volume.imSz1];
handles.display.images.yLimits = [1 handles.volume.imSz2];
handles.display.images.zLimits = [1 handles.volume.imSz3];
handles.display.images.xlimPad = [-0.5 0.5];
handles.display.images.ylimPad = [-0.5 0.5];
handles.display.images.zoomFctr = 1;

set(handles.axes.a1,'xlim',[1 handles.volume.imSz1]+handles.display.images.xlimPad,...
    'ylim',[1 handles.volume.imSz2]+handles.display.images.ylimPad,'Layer','top')
set(handles.axes.a2,'xlim',[1 handles.volume.imSz1]+handles.display.images.xlimPad,...
    'ylim',[1 handles.volume.imSz3]+handles.display.images.ylimPad,'Layer','top')
set(handles.axes.a3,'xlim',[1 handles.volume.imSz2]+handles.display.images.xlimPad,...
    'ylim',[1 handles.volume.imSz3]+handles.display.images.ylimPad,'Layer','top')

%INITIALIZE SYNCED CROSSHAIRS
initializeCrosshair()
%create colored borders around axes
set(handles.axes.a1,'XColor',handles.display.crosshairs.crossCol3,...
    'YColor',handles.display.crosshairs.crossCol3);
set(handles.axes.a2,'XColor',handles.display.crosshairs.crossCol1,...
    'YColor',handles.display.crosshairs.crossCol1);
set(handles.axes.a3,'XColor',handles.display.crosshairs.crossCol2,...
    'YColor',handles.display.crosshairs.crossCol2);

%INITIALIZE SLICE ARROW MARKERS
initializeArrows()

%if IOP defined and volume is aligned to some scanner axes, label orientation in figure
if isfield(handles.volume,'IOP')
    labelPatientOrientation()
else
    delete(findobj('Tag','IOPt1'))
    delete(findobj('Tag','IOPt2'))
    delete(findobj('Tag','IOPt3'))
    delete(findobj('Tag','IOPt4'))
    delete(findobj('Tag','IOPt5'))
    delete(findobj('Tag','IOPt6'))
end
a1T = get(handles.axes.a1,'Title');
a2T = get(handles.axes.a2,'Title');
a3T = get(handles.axes.a3,'Title');
fs = 0.01 * figSz(1) * handles.figure.DPIscl;
if handles.figure.DPIscl~=1
    fs = fs * 0.7;
end
set(a1T,'String',handles.axes.a1title,'FontSize',fs)
set(a2T,'String',handles.axes.a2title,'FontSize',fs)
set(a3T,'String',handles.axes.a3title,'FontSize',fs)

%DEFINE DISPLAY WINDOW AND SLIDERS
if ~handles.display.images.fusedImages
    if get(handles.display.sliders.s4,'Value')==0 && get(handles.display.sliders.s5,'Value')==0
        cLimMin = floor(double(cLimMin));
        cLimMax = ceil(double(cLimMax));
        cLim = [cLimMin cLimMax];
        if ~strcmpi(handles.volume.modality,'PET')
            setCM(gray)
            if strcmpi(handles.volume.modality,'CT')
                cLim = getCTwindow('Body'); %start with body window
                temp = getCTwindow('Bone');
                if cLimMax<temp(2)
                    cLimMax = temp(2);
                end
                temp = getCTwindow('Lung');
                if cLimMin>temp(1)
                    cLimMin = temp(1);
                end
            end
        end
    else
        cLim = [get(handles.display.sliders.s5,'Value') get(handles.display.sliders.s4,'Value')];
        if cLim(1)<cLimMin, cLim(1) = cLimMin; end
        if cLim(1)>cLimMax, cLim(1) = cLimMin; end
        if cLim(2)<cLimMin, cLim(2) = cLimMax; end
        if cLim(2)>cLimMax, cLim(2) = cLimMax; end
    end
    
    int = 100;
    
    set(handles.display.sliders.s4,'Min',cLimMin+diff([cLimMin cLimMax])/int,'Max',cLimMax)
    set(handles.display.sliders.s5,'Min',cLimMin,'Max',cLimMax-diff([cLimMin cLimMax])/int)
    set(handles.display.sliders.s4,'Value',cLim(2))
    set(handles.display.sliders.s5,'Value',cLim(1))
    
    %     %change slider step to accomodate only round numbers
    %     temp1 = get(handles.display.sliders.s4,'SliderStep');
    %     temp2 = ceil(cLimMax)-get(handles.display.sliders.s4,'Min');
    %     if temp1(1)*temp2<1
    %         set(handles.display.sliders.s4,'SliderStep',[1/temp2 10/temp2])
    %         set(handles.display.sliders.s5,'SliderStep',[1/temp2 10/temp2])
    %     end
    
    %pre-calculate slider value values, at 1% steps
    int = 100;
    handles.display.sliders.s4steps = single(...
        (get(handles.display.sliders.s4,'Min'):...
        (get(handles.display.sliders.s4,'Max')-get(handles.display.sliders.s4,'Min'))/(int-1):...
        get(handles.display.sliders.s4,'Max')));
    handles.display.sliders.s5steps = single(...
        (get(handles.display.sliders.s5,'Min'):...
        (get(handles.display.sliders.s5,'Max')-get(handles.display.sliders.s5,'Min'))/(int-1):...
        get(handles.display.sliders.s5,'Max')));
    if get(handles.display.sliders.s4,'Max')-get(handles.display.sliders.s4,'Min')>int
        handles.display.sliders.s4steps = round(handles.display.sliders.s4steps);
    end
    if get(handles.display.sliders.s5,'Max')-get(handles.display.sliders.s5,'Min')>int
        handles.display.sliders.s5steps = round(handles.display.sliders.s5steps);
    end
    
    %text formatting
    s4tString = formatValueString(cLim(2)/handles.display.images.sclFctr);
    s5tString = formatValueString(cLim(1)/handles.display.images.sclFctr);
    
    set(handles.display.sliders.s4t,'String',s4tString)
    set(handles.display.sliders.s5t,'String',s5tString)
    handles.display.images.cLim = cLim;
    
    setDisplayWindow()
else
    if get(handles.display.sliders.s4_1,'Value')==0 && get(handles.display.sliders.s5_1,'Value')==0
        %initiate viewing window on first load
        cLimMin_1 = floor(double(cLimMin_1));
        cLimMax_1 = ceil(double(cLimMax_1));
        cLim_1 = [cLimMin_1 cLimMax_1];
    else
        cLim_1 = [get(handles.display.sliders.s5_1,'Value') get(handles.display.sliders.s4_1,'Value')];
        if cLim_1(1)<cLimMin_1, cLim_1(1) = cLimMin_1; end
        if cLim_1(1)>cLimMax_1, cLim_1(1) = cLimMin_1; end
        if cLim_1(2)<cLimMin_1, cLim_1(2) = cLimMax_1; end
        if cLim_1(2)>cLimMax_1, cLim_1(2) = cLimMax_1; end
    end
    cLimMax_1 = cLim_1(2);
    cLimMin_1 = cLim_1(1);
    
    %always initiate CT with body window
    if strcmpi(handles.volume.modality_1,'CT')
        cLim_1 = getCTwindow('Body'); %start with body window
        temp = getCTwindow('Bone');
        if cLimMax_1<temp(2)
            cLimMax_1 = temp(2);
        end
        temp = getCTwindow('Lung');
        if cLimMin_1>temp(1)
            cLimMin_1 = temp(1);
        end
    end
    
    if get(handles.display.sliders.s4_2,'Value')==0 && get(handles.display.sliders.s5_2,'Value')==0
        %initiate viewing window on first load
        cLimMin_2 = floor(double(cLimMin_2));
        cLimMax_2 = ceil(double(cLimMax_2));
        cLim_2 = [cLimMin_2 cLimMax_2];
    else
        cLim_2 = [get(handles.display.sliders.s5_2,'Value') get(handles.display.sliders.s4_2,'Value')];
        if cLim_2(1)<cLimMin_2, cLim_2(1) = floor(cLimMin_2); end
        if cLim_2(1)>cLimMax_2, cLim_2(1) = floor(cLimMin_2); end
        if cLim_2(2)<cLimMin_2, cLim_2(2) = ceil(cLimMax_2); end
        if cLim_2(2)>cLimMax_2, cLim_2(2) = ceil(cLimMax_2); end
    end
    cLimMax_2 = cLim_2(2);
    cLimMin_2 = cLim_2(1);
    
    %always initiate CT with body window
    if strcmpi(handles.volume.modality_2,'CT')
        cLim_2 = getCTwindow('Body'); %start with body window
        temp = getCTwindow('Bone');
        if cLimMax_2<temp(2)
            cLimMax_2 = temp(2);
        end
        temp = getCTwindow('Lung');
        if cLimMin_2>temp(1)
            cLimMin_2 = temp(1);
        end
    end
    
    int = 100;
    
    set(handles.display.sliders.s4_1,'Min',cLimMin_1+diff([cLimMin_1 cLimMax_1])/int,'Max',cLimMax_1)
    set(handles.display.sliders.s5_1,'Min',cLimMin_1,'Max',cLimMax_1-diff([cLimMin_1 cLimMax_1])/int)
    set(handles.display.sliders.s4_2,'Min',cLimMin_2+diff([cLimMin_2 cLimMax_2])/int,'Max',cLimMax_2)
    set(handles.display.sliders.s5_2,'Min',cLimMin_2,'Max',cLimMax_2-diff([cLimMin_2 cLimMax_2])/int)
    set(handles.display.sliders.s4_1,'Value',cLim_1(2))
    set(handles.display.sliders.s5_1,'Value',cLim_1(1))
    set(handles.display.sliders.s4_2,'Value',cLim_2(2))
    set(handles.display.sliders.s5_2,'Value',cLim_2(1))
    %     set(handles.display.sliders.s4_1,'Value',round(cLim_1(2)))
    %     set(handles.display.sliders.s5_1,'Value',round(cLim_1(1)))
    %     set(handles.display.sliders.s4_2,'Value',round(cLim_2(2)))
    %     set(handles.display.sliders.s5_2,'Value',round(cLim_2(1)))
    %     %change slider steps to accomodate only round numbers
    %     temp1 = get(handles.display.sliders.s4_1,'SliderStep');
    %     temp2 = get(handles.display.sliders.s4_1,'Max')-get(handles.display.sliders.s4_1,'Min');
    %     if temp1(1)*temp2<1
    %         set(handles.display.sliders.s4_1,'SliderStep',[1/temp2 10/temp2])
    %         set(handles.display.sliders.s5_1,'SliderStep',[1/temp2 10/temp2])
    %     end
    %     temp1 = get(handles.display.sliders.s4_2,'SliderStep');
    %     temp2 = get(handles.display.sliders.s4_2,'Max')-get(handles.display.sliders.s4_2,'Min');
    %     if temp1(1)*temp2<1
    %         set(handles.display.sliders.s4_2,'SliderStep',[1/temp2 10/temp2])
    %         set(handles.display.sliders.s5_2,'SliderStep',[1/temp2 10/temp2])
    %     end
    
    %pre-calculate slider value values, at 1% steps
    int = 100;
    handles.display.sliders.s4steps_1 = single(...
        (get(handles.display.sliders.s4_1,'Min'):...
        (get(handles.display.sliders.s4_1,'Max')-get(handles.display.sliders.s4_1,'Min'))/(int-1):...
        get(handles.display.sliders.s4_1,'Max')));
    handles.display.sliders.s5steps_1 = single(...
        (get(handles.display.sliders.s5_1,'Min'):...
        (get(handles.display.sliders.s5_1,'Max')-get(handles.display.sliders.s5_1,'Min'))/(int-1):...
        get(handles.display.sliders.s5_1,'Max')));
    handles.display.sliders.s4steps_2 = single(...
        (get(handles.display.sliders.s4_2,'Min'):...
        (get(handles.display.sliders.s4_2,'Max')-get(handles.display.sliders.s4_2,'Min'))/(int-1):...
        get(handles.display.sliders.s4_2,'Max')));
    handles.display.sliders.s5steps_2 = single(...
        (get(handles.display.sliders.s5_2,'Min'):...
        (get(handles.display.sliders.s5_2,'Max')-get(handles.display.sliders.s5_2,'Min'))/(int-1):...
        get(handles.display.sliders.s5_2,'Max')));
    if get(handles.display.sliders.s4_1,'Max')-get(handles.display.sliders.s4_1,'Min')>int
        handles.display.sliders.s4steps_1 = round(handles.display.sliders.s4steps_1);
    end
    if get(handles.display.sliders.s5_1,'Max')-get(handles.display.sliders.s5_1,'Min')>int
        handles.display.sliders.s5steps_1 = round(handles.display.sliders.s5steps_1);
    end
    if get(handles.display.sliders.s4_2,'Max')-get(handles.display.sliders.s4_2,'Min')>int
        handles.display.sliders.s4steps_2 = round(handles.display.sliders.s4steps_2);
    end
    if get(handles.display.sliders.s5_2,'Max')-get(handles.display.sliders.s5_2,'Min')>int
        handles.display.sliders.s5steps_2 = round(handles.display.sliders.s5steps_2);
    end
    
    %text formatting
    s4tString_1 = formatValueString(cLim_1(2)/handles.display.images.sclFctr_1);
    s5tString_1 = formatValueString(cLim_1(1)/handles.display.images.sclFctr_1);
    s4tString_2 = formatValueString(cLim_2(2)/handles.display.images.sclFctr_2);
    s5tString_2 = formatValueString(cLim_2(1)/handles.display.images.sclFctr_2);
    
    set(handles.display.sliders.s4t_1,'String',s4tString_1)
    set(handles.display.sliders.s5t_1,'String',s5tString_1)
    set(handles.display.sliders.s4t_2,'String',s4tString_2)
    set(handles.display.sliders.s5t_2,'String',s5tString_2)
    
    handles.display.images.cLim_1 = cLim_1;
    handles.display.images.cLim_2 = cLim_2;
    
    setFusedDisplayWindow()
end

%INITIATE SLICE NUMBER TEXT BOXES BELOW SLIDER BAR HANDLE
s1Pos = get(handles.display.sliders.s1,'Position');
s2Pos = get(handles.display.sliders.s2,'Position');
s3Pos = get(handles.display.sliders.s3,'Position');
handles.display.sliders.stSz = [0.03*figSz(1) handles.display.sliders.sldrWidth];
handles.display.sliders.st_y = s1Pos(2)-1.3*handles.display.sliders.sldrWidth; %ASSUMES SAME SLIDER Y POSITION
arrowButtonWidth = 16; %TO ACCOUNT FOR SLIDER "TROUGH" SIZE
%DEFINE NAVIGATION SLIDERS AND COMPUTE ARRRAYS OF TEXT BOX LOCATIONS
if handles.volume.imSz3>1
    set(handles.display.sliders.s1,'Min',1,'Max',handles.volume.imSz3,'Value',handles.display.images.z,...
        'SliderStep',[1/(handles.volume.imSz3-1) 5/(handles.volume.imSz3-1)])
    s1hw = getSliderKnobWidth(ax1Pos(3),arrowButtonWidth,get(handles.display.sliders.s1,'SliderStep'));
    s1m = (s1Pos(3)-s1hw-2*arrowButtonWidth)/(handles.volume.imSz3-1);
    s1b = s1Pos(1)+arrowButtonWidth+s1hw/2-...
        ((s1Pos(3)-2*arrowButtonWidth-s1hw)/(handles.volume.imSz3-1));
    handles.display.sliders.s1tLocs = ...
        s1m.*(1:handles.volume.imSz3)+s1b-handles.display.sliders.stSz(1)/2;
else
    set(handles.display.sliders.s1,'Min',1,'Max',handles.volume.imSz3,'Value',handles.display.images.z)
    handles.display.sliders.s1tLocs = s1Pos(1)+s1Pos(3)/2-handles.display.sliders.stSz(1)/2;
end
if handles.volume.imSz2>1
    set(handles.display.sliders.s2,'Min',1,'Max',handles.volume.imSz2,'Value',handles.display.images.y,...
        'SliderStep',[1/(handles.volume.imSz2-1) 5/(handles.volume.imSz2-1)])
    s2hw = getSliderKnobWidth(ax2Pos(3),arrowButtonWidth,get(handles.display.sliders.s2,'SliderStep'));
    s2m = (s2Pos(3)-s2hw-2*arrowButtonWidth)/(handles.volume.imSz2-1);
    s2b = s2Pos(1)+arrowButtonWidth+s2hw/2-...
        ((s2Pos(3)-2*arrowButtonWidth-s2hw)/(handles.volume.imSz2-1));
    handles.display.sliders.s2tLocs = ...
        s2m.*(1:handles.volume.imSz2)+s2b-handles.display.sliders.stSz(1)/2;
else
    set(handles.display.sliders.s2,'Min',1,'Max',handles.volume.imSz2,'Value',handles.display.images.y)
    handles.display.sliders.s2tLocs = s2Pos(1)+s2Pos(3)/2-handles.display.sliders.stSz(1)/2;
end
if handles.volume.imSz1>1
    set(handles.display.sliders.s3,'Min',1,'Max',handles.volume.imSz1,'Value',handles.display.images.x,...
        'SliderStep',[1/(handles.volume.imSz1-1) 5/(handles.volume.imSz1-1)])
    s3hw = getSliderKnobWidth(ax3Pos(3),arrowButtonWidth,get(handles.display.sliders.s3,'SliderStep'));
    s3m = (s3Pos(3)-s3hw-2*arrowButtonWidth)/(handles.volume.imSz1-1);
    s3b = s3Pos(1)+arrowButtonWidth+s3hw/2-...
        ((s3Pos(3)-2*arrowButtonWidth-s3hw)/(handles.volume.imSz1-1));
    handles.display.sliders.s3tLocs = ...
        s3m.*(1:handles.volume.imSz1)+s3b-handles.display.sliders.stSz(1)/2;
else
    set(handles.display.sliders.s3,'Min',1,'Max',handles.volume.imSz1,'Value',handles.display.images.x)
    handles.display.sliders.s3tLocs = s3Pos(1)+s3Pos(3)/2-handles.display.sliders.stSz(1)/2;
end
%NOW DO TEXT BOXES
s1tPos = [handles.display.sliders.s1tLocs(get(handles.display.sliders.s1,'Value')) ...
    handles.display.sliders.st_y handles.display.sliders.stSz];
s2tPos = [handles.display.sliders.s2tLocs(get(handles.display.sliders.s2,'Value')) ...
    handles.display.sliders.st_y handles.display.sliders.stSz];
s3tPos = [handles.display.sliders.s3tLocs(get(handles.display.sliders.s3,'Value')) ...
    handles.display.sliders.st_y handles.display.sliders.stSz];
stFS = 0.007*figSz(1)*handles.figure.DPIscl;

str = ['''Style'',''Edit'','...
    '''Units'',''Pixel'','...
    '''Position'',[' num2str(s1tPos) '],'...
    '''FontSize'',' num2str(stFS) ','...
    '''HorizontalAlignment'',''Center'','...
    '''FontWeight'',''Bold'','...
    '''Visible'',''Off'','...
    '''String'',''' num2str(handles.display.images.z) ''','...
    '''Callback'',@a1TextCallback,'...
    '''BackgroundColor'',[' num2str(get(findobj('Tag','3DviewerFig'),'color')) ']'];
handles.display.sliders.s1t = initiateObject('s1Text','uicontrol',str);

str = ['''Style'',''Edit'','...
    '''Units'',''Pixel'','...
    '''Position'',[' num2str(s2tPos) '],'...
    '''FontSize'',' num2str(stFS) ','...
    '''HorizontalAlignment'',''Center'','...
    '''FontWeight'',''Bold'','...
    '''Visible'',''Off'','...
    '''String'',''' num2str(handles.display.images.y) ''','...
    '''Callback'',@a2TextCallback,'...
    '''BackgroundColor'',[' num2str(get(findobj('Tag','3DviewerFig'),'color')) ']'];
handles.display.sliders.s2t = initiateObject('s2Text','uicontrol',str);

str = ['''Style'',''Edit'','...
    '''Units'',''Pixel'','...
    '''Position'',[' num2str(s3tPos) '],'...
    '''FontSize'',' num2str(stFS) ','...
    '''HorizontalAlignment'',''Center'','...
    '''FontWeight'',''Bold'','...
    '''Visible'',''Off'','...
    '''String'',''' num2str(handles.display.images.x) ''','...
    '''Callback'',@a3TextCallback,'...
    '''BackgroundColor'',[' num2str(get(findobj('Tag','3DviewerFig'),'color')) ']'];
handles.display.sliders.s3t = initiateObject('s3Text','uicontrol',str);

%SET BOOLEAN FOR CURSOR OVER IMAGES
handles.display.cursor.overAxes = false;

%STAMP FIGURE WITH IMAGE VOLUME INFORMATION
stampImageInfo()

drawnow

end

function h = setVolume(h,a)

delete(findobj('Tag','a1im_1'))
delete(findobj('Tag','a1im_2'))
delete(findobj('Tag','a2im_1'))
delete(findobj('Tag','a2im_2'))
delete(findobj('Tag','a3im_1'))
delete(findobj('Tag','a3im_2'))

global handles

switch a
    case handles.axes.a1
        if ~isempty(h)
            showImage1(handles.volume.vol,h,handles.display.images.z,handles.display.images.fltr)
        else
            axes(a)
            h = imagesc(rot90(squeeze(handles.volume.vol(:,:,handles.display.images.z))),...
                'HitTest','off','Visible','off','Tag','a1im');
%             h = imagesc(imfilter(rot90(squeeze(handles.volume.vol(:,:,handles.display.images.z))),handles.display.images.fltr,'replicate'),...
%                 'HitTest','off','Visible','off','Tag','a1im');
            hold on
        end
    case handles.axes.a2
        if ~isempty(h)
            showImage2(handles.volume.vol,h,handles.display.images.y,handles.display.images.fltr)
        else
            axes(a)
            h = imagesc(imfilter(permute(squeeze(handles.volume.vol(:,handles.display.images.y,:)),[2 1]),handles.display.images.fltr,'replicate'),...
                'HitTest','off','Visible','off','Tag','a2im');
%             h = imagesc(imfilter(flipud(rot90(squeeze(handles.volume.vol(:,handles.display.images.y,:)))),handles.display.images.fltr,'replicate'),...
%                 'HitTest','off','Visible','off','Tag','a2im');
            hold on
        end
    case handles.axes.a3
        if ~isempty(h)
            showImage3(handles.volume.vol,h,handles.display.images.x,handles.display.images.fltr)
        else
            axes(a)
            h = imagesc(imfilter(permute(squeeze(handles.volume.vol(handles.display.images.x,:,:)),[2 1]),handles.display.images.fltr,'replicate'),...
                'HitTest','off','Visible','off','Tag','a3im');
%             h = imagesc(imfilter(flipud(rot90(squeeze(handles.volume.vol(handles.display.images.x,:,:)))),handles.display.images.fltr,'replicate'),...
%                 'HitTest','off','Visible','off','Tag','a3im');
            hold on
        end
end

end

function [h1, h2] = setFusedVolume(h1,h2,a,clim1,clim2,cm1,cm2)

delete(findobj('Tag','a1im'))
delete(findobj('Tag','a2im'))
delete(findobj('Tag','a3im'))

%h1 is top image

global handles

switch a
    case handles.axes.a1
        if ~isempty(h1) && ~isempty(h2)
            showFusedImage1(a,handles.volume.vol_1,handles.volume.vol_2,h1,h2,handles.display.images.z,clim1,clim2,...
                handles.display.images.cm_1,handles.display.images.cm_2,handles.display.images.fltr)            
            alpha(h1,1-get(handles.display.sliders.s6_1,'Value'));
            alpha(h2,1);
        else
            axes(a)
            h2 = imagesc(rot90(squeeze(handles.volume.vol_2(:,:,handles.display.images.z))),...
                'HitTest','off','Tag','a1im_2');
%             h2 = imagesc(imfilter(rot90(squeeze(handles.volume.vol_2(:,:,handles.display.images.z))),handles.display.images.fltr,'replicate'),...
%                 'HitTest','off','Tag','a1im_2');
            colormap(cm2)
            set(a,'CLim',clim2)
            freezeImageColors(a,h2)
            h1 = imagesc(rot90(squeeze(handles.volume.vol_1(:,:,handles.display.images.z))),...
                'HitTest','off','Tag','a1im_1');
%             h1 = imagesc(imfilter(rot90(squeeze(handles.volume.vol_1(:,:,handles.display.images.z))),handles.display.images.fltr,'replicate'),...
%                 'HitTest','off','Tag','a1im_1');
            colormap(cm1)
            set(a,'CLim',clim1)
            
            set(h1,'Visible','off')
            set(h2,'Visible','off')
            
            alpha(h1,1-get(handles.display.sliders.s6_1,'Value'));
            alpha(h2,1);
            
        end
    case handles.axes.a2
        if ~isempty(h1) && ~isempty(h2)
            showFusedImage2(a,handles.volume.vol_1,handles.volume.vol_2,h1,h2,handles.display.images.y,clim1,clim2,...
                handles.display.images.cm_1,handles.display.images.cm_2,handles.display.images.fltr)            
            alpha(h1,1-get(handles.display.sliders.s6_1,'Value'));
            alpha(h2,1);
        else
            axes(a)            
            h2 = imagesc(permute(squeeze(handles.volume.vol_2(:,handles.display.images.y,:)),[2 1]),...
                'HitTest','off','Tag','a2im_2');
%             h2 = imagesc(imfilter(flipud(rot90(squeeze(handles.volume.vol_2(:,handles.display.images.y,:)))),handles.display.images.fltr,'replicate'),...
%                 'HitTest','off','Tag','a2im_2');
            colormap(cm2)
            set(a,'CLim',clim2)
            freezeImageColors(a,h2)
            h1 = imagesc(permute(squeeze(handles.volume.vol_1(:,handles.display.images.y,:)),[2 1]),...
                'HitTest','off','Tag','a2im_1');
%             h1 = imagesc(imfilter(flipud(rot90(squeeze(handles.volume.vol_1(:,handles.display.images.y,:)))),handles.display.images.fltr,'replicate'),...
%                 'HitTest','off','Tag','a2im_1');
            colormap(cm1)
            set(a,'CLim',clim1)
            
            set(h1,'Visible','off')
            set(h2,'Visible','off')
            
            alpha(h1,1-get(handles.display.sliders.s6_1,'Value'));
            alpha(h2,1);
            
        end
    case handles.axes.a3
        if ~isempty(h1) && ~isempty(h2)
            showFusedImage3(a,handles.volume.vol_1,handles.volume.vol_2,h1,h2,handles.display.images.x,clim1,clim2,...
                handles.display.images.cm_1,handles.display.images.cm_2,handles.display.images.fltr)            
            alpha(h1,1-get(handles.display.sliders.s6_1,'Value'));
            alpha(h2,1);
        else
            axes(a)
            h2 = imagesc(permute(squeeze(handles.volume.vol_2(handles.display.images.x,:,:)),[2 1]),...
                'HitTest','off','Tag','a3im_2');
%             h2 = imagesc(imfilter(flipud(rot90(squeeze(handles.volume.vol_2(handles.display.images.x,:,:)))),handles.display.images.fltr,'replicate'),...
%                 'HitTest','off','Tag','a3im_2');
            colormap(cm2)
            set(a,'CLim',clim2)
            freezeImageColors(a,h2)
            h1 = imagesc(permute(squeeze(handles.volume.vol_1(handles.display.images.x,:,:)),[2 1]),...
                'HitTest','off','Tag','a3im_1');
%             h1 = imagesc(imfilter(flipud(rot90(squeeze(handles.volume.vol_1(handles.display.images.x,:,:)))),handles.display.images.fltr,'replicate'),...
%                 'HitTest','off','Tag','a3im_1');
            colormap(cm1)
            set(a,'CLim',clim1)
            
            set(h1,'Visible','off')
            set(h2,'Visible','off')
            
            alpha(h1,1-get(handles.display.sliders.s6_1,'Value'));
            alpha(h2,1);
            
        end
end

end

function showImage1(vol,h,idx,fltr)

set(h,'CData',rot90(squeeze(vol(:,:,idx))))
% set(h,'CData',imfilter(rot90(squeeze(vol(:,:,idx))),fltr,'replicate'))

end

function showImage2(vol,h,idx,fltr)

set(h,'CData',permute(squeeze(vol(:,idx,:)),[2 1]))
% set(h,'CData',imfilter(permute(squeeze(vol(:,idx,:)),[2 1]),fltr,'replicate'))

end

function showImage3(vol,h,idx,fltr)

set(h,'CData',permute(squeeze(vol(idx,:,:)),[2 1]))
% set(h,'CData',imfilter(permute(squeeze(vol(idx,:,:)),[2 1]),fltr,'replicate'))

end

function showFusedImage1(a,vol1,vol2,h1,h2,idx,clim1,clim2,cm1,cm2,fltr)

colormap(cm2)
set(a,'CLim',clim2)
set(h2,'cdata',rot90(squeeze(vol2(:,:,idx))))
% set(h2,'cdata',imfilter(rot90(squeeze(vol2(:,:,idx))),fltr,'replicate'))

freezeImageColors(a,h2)
colormap(cm1)
set(a,'CLim',clim1)
set(h1,'cdata',rot90(squeeze(vol1(:,:,idx))))
% set(h1,'cdata',imfilter(rot90(squeeze(vol1(:,:,idx))),fltr,'replicate'))

end

function showFusedImage2(a,vol1,vol2,h1,h2,idx,clim1,clim2,cm1,cm2,fltr)

colormap(cm2)
set(a,'CLim',clim2)
set(h2,'cdata',permute(squeeze(vol2(:,idx,:)),[2 1]))
% set(h2,'cdata',imfilter(permute(squeeze(vol2(:,idx,:)),[2 1]),fltr,'replicate'))

freezeImageColors(a,h2)
colormap(cm1)
set(a,'CLim',clim1)
set(h1,'cdata',permute(squeeze(vol1(:,idx,:)),[2 1]))
% set(h1,'cdata',imfilter(permute(squeeze(vol1(:,idx,:)),[2 1]),fltr,'replicate'))

end

function showFusedImage3(a,vol1,vol2,h1,h2,idx,clim1,clim2,cm1,cm2,fltr)

colormap(cm2)
set(a,'CLim',clim2)
set(h2,'cdata',permute(squeeze(vol2(idx,:,:)),[2 1]))
% set(h2,'cdata',imfilter(permute(squeeze(vol2(idx,:,:)),[2 1]),fltr,'replicate'))

freezeImageColors(a,h2)
colormap(cm1)
set(a,'CLim',clim1)
set(h1,'cdata',permute(squeeze(vol1(idx,:,:)),[2 1]))
% set(h1,'cdata',imfilter(permute(squeeze(vol1(idx,:,:)),[2 1]),fltr,'replicate'))

end

function setDisplaySliders

global handles

if get(handles.display.sliders.s4,'Value')<=handles.display.images.cLim(1)
    set(handles.display.sliders.s4,'Value',handles.display.images.cLim(1)+1)
    handles.display.images.cLim(2) = get(handles.display.sliders.s4,'Value');
end
if get(handles.display.sliders.s5,'Value')>=handles.display.images.cLim(2)
    set(handles.display.sliders.s5,'Value',handles.display.images.cLim(2)-1)
    handles.display.images.cLim(1) = get(handles.display.sliders.s5,'Value');
end

end

function setDisplaySliders_1

global handles

if get(handles.display.sliders.s4_1,'Value')<=handles.display.images.cLim_1(1)
    set(handles.display.sliders.s4_1,'Value',handles.display.images.cLim_1(1)+1)
    handles.display.images.cLim_1(2) = get(handles.display.sliders.s4_1,'Value');
end
if get(handles.display.sliders.s5_1,'Value')>=handles.display.images.cLim_1(2)
    set(handles.display.sliders.s5_1,'Value',handles.display.images.cLim_1(2)-1)
    handles.display.images.cLim_1(1) = get(handles.display.sliders.s5_1,'Value');
end

end

function setDisplaySliders_2

global handles

if get(handles.display.sliders.s4_2,'Value')<=handles.display.images.cLim_2(1)
    set(handles.display.sliders.s4_2,'Value',handles.display.images.cLim_2(1)+1)
    handles.display.images.cLim_2(2) = get(handles.display.sliders.s4_2,'Value');
end
if get(handles.display.sliders.s5_2,'Value')>=handles.display.images.cLim_2(2)
    set(handles.display.sliders.s5_2,'Value',handles.display.images.cLim_2(2)-1)
    handles.display.images.cLim_2(1) = get(handles.display.sliders.s5_2,'Value');
end

end

function setDisplayWindow

global handles

set(handles.axes.a1,'Clim',handles.display.images.cLim)
set(handles.axes.a2,'Clim',handles.display.images.cLim)
set(handles.axes.a3,'Clim',handles.display.images.cLim)

end

function setFusedDisplayWindow

global handles

showFusedImage1(handles.axes.a1,handles.volume.vol_1,handles.volume.vol_2,...
    handles.display.images.h1_1,handles.display.images.h1_2,handles.display.images.z,...
    handles.display.images.cLim_1,handles.display.images.cLim_2,...
    handles.display.images.cm_1,handles.display.images.cm_2,handles.display.images.fltr)
showFusedImage2(handles.axes.a2,handles.volume.vol_1,handles.volume.vol_2,...
    handles.display.images.h2_1,handles.display.images.h2_2,handles.display.images.y,...
    handles.display.images.cLim_1,handles.display.images.cLim_2,...
    handles.display.images.cm_1,handles.display.images.cm_2,handles.display.images.fltr)
showFusedImage3(handles.axes.a3,handles.volume.vol_1,handles.volume.vol_2,...
    handles.display.images.h3_1,handles.display.images.h3_2,handles.display.images.x,...
    handles.display.images.cLim_1,handles.display.images.cLim_2,...
    handles.display.images.cm_1,handles.display.images.cm_2,handles.display.images.fltr)

end

function initializeCrosshair

global handles

Xlim = xlim(handles.axes.a1);
Ylim = ylim(handles.axes.a1);
Zlim = ylim(handles.axes.a2);

handles.display.crosshairs.crossGap = 15;  %mm
% handles.display.crosshairs.crossCol1 = 'b';
% handles.display.crosshairs.crossCol2 = 'g';
% handles.display.crosshairs.crossCol3 = 'r';
handles.display.crosshairs.crossCol1 = [0 0.4470 0.7411];
handles.display.crosshairs.crossCol2 = [0.4460 0.6740 0.1880];
handles.display.crosshairs.crossCol3 = [0.8500 0.3250 0.0980];

test = findobj('Tag','a1crossH1');

if ~isempty(test) %ASSUME IF ONE FOUND, ALL EXIST
    %a1
    crossX = handles.display.images.x;
    crossY = handles.display.images.y;
    xGap = round(handles.display.crosshairs.crossGap/handles.volume.voxDimX/2);
    yGap = round(handles.display.crosshairs.crossGap/handles.volume.voxDimY/2);
    handles.display.crosshairs.a1crossH1 = test;
    set(handles.display.crosshairs.a1crossH1,'XData',(Xlim(1)-1:crossX-xGap),'YData',crossY*(ones(size((Xlim(1)-1:crossX-xGap)))),'Visible','off')
    handles.display.crosshairs.a1crossH2 = findobj('Tag','a1crossH2');
    set(handles.display.crosshairs.a1crossH2,'XData',(crossX+xGap:Xlim(2)+1),'YData',crossY*(ones(size((crossX+xGap:Xlim(2)+1)))),'Visible','off')
    handles.display.crosshairs.a1crossV1 = findobj('Tag','a1crossV1');
    set(handles.display.crosshairs.a1crossV1,'XData',crossX*(ones(size((Ylim(1)-1:crossY-yGap)))),'YData',(Ylim(1)-1:crossY-yGap),'Visible','off')
    handles.display.crosshairs.a1crossV2 = findobj('Tag','a1crossV2');
    set(handles.display.crosshairs.a1crossV2,'XData',crossX*(ones(size((crossY+yGap:Ylim(2)+1)))),'YData',(crossY+yGap:Ylim(2)+1),'Visible','off')
    %a2
    crossY = handles.display.images.z;
    xGap = round(handles.display.crosshairs.crossGap/handles.volume.voxDimX/2);
    yGap = round(handles.display.crosshairs.crossGap/handles.volume.voxDimZ/2);
    handles.display.crosshairs.a2crossH1 = findobj('Tag','a2crossH1');
    set(handles.display.crosshairs.a2crossH1,'XData',(Xlim(1)-1:crossX-xGap),'YData',crossY*(ones(size((Xlim(1)-1:crossX-xGap)))),'Visible','off')
    handles.display.crosshairs.a2crossH2 = findobj('Tag','a2crossH2');
    set(handles.display.crosshairs.a2crossH2,'XData',(crossX+xGap:Xlim(2)+1),'YData',crossY*(ones(size((crossX+xGap:Xlim(2)+1)))),'Visible','off')
    handles.display.crosshairs.a2crossV1 = findobj('Tag','a2crossV1');
    set(handles.display.crosshairs.a2crossV1,'XData',crossX*(ones(size((Zlim(1)-1:crossY-yGap)))),'YData',(Zlim(1)-1:crossY-yGap),'Visible','off')
    handles.display.crosshairs.a2crossV2 = findobj('Tag','a2crossV2');
    set(handles.display.crosshairs.a2crossV2,'XData',crossX*(ones(size((crossY+yGap:Zlim(2)+1)))),'YData',(crossY+yGap:Zlim(2)+1),'Visible','off')
    %a3
    crossX = handles.display.images.y;
    xGap = round(handles.display.crosshairs.crossGap/handles.volume.voxDimY/2);
    yGap = round(handles.display.crosshairs.crossGap/handles.volume.voxDimZ/2);
    handles.display.crosshairs.a3crossH1 = findobj('Tag','a3crossH1');
    set(handles.display.crosshairs.a3crossH1,'XData',(Ylim(1)-1:crossX-xGap),'YData',crossY*(ones(size((Ylim(1)-1:crossX-xGap)))),'Visible','off')
    handles.display.crosshairs.a3crossH2 = findobj('Tag','a3crossH2');
    set(handles.display.crosshairs.a3crossH2,'XData',(crossX+xGap:Ylim(2)+1),'YData',crossY*(ones(size((crossX+xGap:Ylim(2)+1)))),'Visible','off')
    handles.display.crosshairs.a3crossV1 = findobj('Tag','a3crossV1');
    set(handles.display.crosshairs.a3crossV1,'XData',crossX*(ones(size((Zlim(1)-1:crossY-yGap)))),'YData',(Zlim(1)-1:crossY-yGap),'Visible','off')
    handles.display.crosshairs.a3crossV2 = findobj('Tag','a3crossV2');
    set(handles.display.crosshairs.a3crossV2,'XData',crossX*(ones(size((crossY+yGap:Zlim(2)+1)))),'YData',(crossY+yGap:Zlim(2)+1),'Visible','off')
else
    %a1
    axes(handles.axes.a1)
    crossX = handles.display.images.x;
    crossY = handles.display.images.y;
    xGap = round(handles.display.crosshairs.crossGap/handles.volume.voxDimX/2);
    yGap = round(handles.display.crosshairs.crossGap/handles.volume.voxDimY/2);
    handles.display.crosshairs.a1crossH1 = plot((Xlim(1)-1:crossX-xGap),crossY*(ones(size((Xlim(1)-1:crossX-xGap)))),...
        'Color',handles.display.crosshairs.crossCol1,'HitTest','off','Visible','off','Tag','a1crossH1');
    handles.display.crosshairs.a1crossH2 = plot((crossX+xGap:Xlim(2)+1),crossY*(ones(size((crossX+xGap:Xlim(2)+1)))),...
        'Color',handles.display.crosshairs.crossCol1,'HitTest','off','Visible','off','Tag','a1crossH2');
    handles.display.crosshairs.a1crossV1 = plot(crossX*(ones(size((Ylim(1)-1:crossY-yGap)))),(Ylim(1)-1:crossY-yGap),...
        'Color',handles.display.crosshairs.crossCol2,'HitTest','off','Visible','off','Tag','a1crossV1');
    handles.display.crosshairs.a1crossV2 = plot(crossX*(ones(size((crossY+yGap:Ylim(2)+1)))),(crossY+yGap:Ylim(2)+1),...
        'Color',handles.display.crosshairs.crossCol2,'HitTest','off','Visible','off','Tag','a1crossV2');
    %a2
    axes(handles.axes.a2)
    crossY = handles.display.images.z;
    xGap = round(handles.display.crosshairs.crossGap/handles.volume.voxDimX/2);
    yGap = round(handles.display.crosshairs.crossGap/handles.volume.voxDimZ/2);
    handles.display.crosshairs.a2crossH1 = plot((Xlim(1)-1:crossX-xGap),crossY*(ones(size((Xlim(1)-1:crossX-xGap)))),...
        'Color',handles.display.crosshairs.crossCol3,'HitTest','off','Visible','off','Tag','a2crossH1');
    handles.display.crosshairs.a2crossH2 = plot((crossX+xGap:Xlim(2)+1),crossY*(ones(size((crossX+xGap:Xlim(2)+1)))),...
        'Color',handles.display.crosshairs.crossCol3,'HitTest','off','Visible','off','Tag','a2crossH2');
    handles.display.crosshairs.a2crossV1 = plot(crossX*(ones(size((Zlim(1)-1:crossY-yGap)))),(Zlim(1)-1:crossY-yGap),...
        'Color',handles.display.crosshairs.crossCol2,'HitTest','off','Visible','off','Tag','a2crossV1');
    handles.display.crosshairs.a2crossV2 = plot(crossX*(ones(size((crossY+yGap:Zlim(2)+1)))),(crossY+yGap:Zlim(2)+1),...
        'Color',handles.display.crosshairs.crossCol2,'HitTest','off','Visible','off','Tag','a2crossV2');
    %a3
    axes(handles.axes.a3)
    crossX = handles.display.images.y;
    xGap = round(handles.display.crosshairs.crossGap/handles.volume.voxDimY/2);
    yGap = round(handles.display.crosshairs.crossGap/handles.volume.voxDimZ/2);
    handles.display.crosshairs.a3crossH1 = plot((Ylim(1)-1:crossX-xGap),crossY*(ones(size((Ylim(1)-1:crossX-xGap)))),...
        'Color',handles.display.crosshairs.crossCol3,'HitTest','off','Visible','off','Tag','a3crossH1');
    handles.display.crosshairs.a3crossH2 = plot((crossX+xGap:Ylim(2)+1),crossY*(ones(size((crossX+xGap:Ylim(2)+1)))),...
        'Color',handles.display.crosshairs.crossCol3,'HitTest','off','Visible','off','Tag','a3crossH2');
    handles.display.crosshairs.a3crossV1 = plot(crossX*(ones(size((Zlim(1)-1:crossY-yGap)))),(Zlim(1)-1:crossY-yGap),...
        'Color',handles.display.crosshairs.crossCol1,'HitTest','off','Visible','off','Tag','a3crossV1');
    handles.display.crosshairs.a3crossV2 = plot(crossX*(ones(size((crossY+yGap:Zlim(2)+1)))),(crossY+yGap:Zlim(2)+1),...
        'Color',handles.display.crosshairs.crossCol1,'HitTest','off','Visible','off','Tag','a3crossV2');
end
handles.display.crosshairs.crossBoo = false;

end

function initializeArrows

global handles

figSz = get(handles.figure.f1,'Position');
figSz = figSz(3:4);

a1Xlim = xlim(handles.axes.a1);
a1Ylim = ylim(handles.axes.a1);
a2Xlim = xlim(handles.axes.a2);
a2Ylim = ylim(handles.axes.a2);
a3Xlim = xlim(handles.axes.a3);
a3Ylim = ylim(handles.axes.a3);

%FIND AXES FIGURE NORMALIZED COORDINATES
units = get(handles.axes.a1,'Units');
set(handles.axes.a1,'Units','pixels')
set(handles.axes.a2,'Units','pixels')
set(handles.axes.a3,'Units','pixels')
figPos = get(handles.figure.f1,'Position');
ax1Pos = get(handles.axes.a1,'Position');
ax2Pos = get(handles.axes.a2,'Position');
ax3Pos = get(handles.axes.a3,'Position');
set(handles.axes.a1,'Units',units)
set(handles.axes.a2,'Units',units)
set(handles.axes.a3,'Units',units)

%normalized axes start and end points
ax1hStart = ax1Pos(1)/figPos(3);
ax1hEnd = (ax1Pos(1)+ax1Pos(3))/figPos(3);
ax1vStart = ax1Pos(2)/figPos(4);
ax1vEnd = (ax1Pos(2)+ax1Pos(4))/figPos(4);
ax2hStart = ax2Pos(1)/figPos(3);
ax2hEnd = (ax2Pos(1)+ax2Pos(3))/figPos(3);
ax2vStart = ax2Pos(2)/figPos(4);
ax2vEnd = (ax2Pos(2)+ax2Pos(4))/figPos(4);
ax3hStart = ax3Pos(1)/figPos(3);
ax3hEnd = (ax3Pos(1)+ax3Pos(3))/figPos(3);
ax3vStart = ax3Pos(2)/figPos(4);
ax3vEnd = (ax3Pos(2)+ax3Pos(4))/figPos(4);

%get normalized size of all axes pixels, accounting for image and padding
voxDim1h = ax1Pos(3)/figPos(3)/diff(a1Xlim);
voxDim1v = ax1Pos(4)/figPos(4)/diff(a1Ylim);
voxDim2h = ax2Pos(3)/figPos(3)/diff(a2Xlim);
voxDim2v = ax2Pos(4)/figPos(4)/diff(a2Ylim);
voxDim3h = ax3Pos(3)/figPos(3)/diff(a3Xlim);
voxDim3v = ax3Pos(4)/figPos(4)/diff(a3Ylim);
%normalized (total, including unplotted) image start and end postions
im1hStart = round((ax1Pos(1)/figPos(3)-((a1Xlim(1)-1)*voxDim1h))*1e8)/1e8;
im1hEnd = round(((ax1Pos(1)+ax1Pos(3))/figPos(3)+((handles.volume.imSz1-a1Xlim(2))*voxDim1h))*1e8)/1e8;
im1vStart = round((ax1Pos(2)/figPos(4)-((a1Ylim(1)-1)*voxDim1v))*1e8)/1e8;
im1vEnd = round(((ax1Pos(2)+ax1Pos(4))/figPos(4)+((handles.volume.imSz2-a1Ylim(2))*voxDim1v))*1e8)/1e8;
im2hStart = round((ax2Pos(1)/figPos(3)-((a2Xlim(1)-1)*voxDim2h))*1e8)/1e8;
im2hEnd = round(((ax2Pos(1)+ax2Pos(3))/figPos(3)+((handles.volume.imSz1-a2Xlim(2))*voxDim2h))*1e8)/1e8;
im2vStart = round((ax2Pos(2)/figPos(4)-((a2Ylim(1)-1)*voxDim2v))*1e8)/1e8;
im2vEnd = round(((ax2Pos(2)+ax2Pos(4))/figPos(4)+((handles.volume.imSz3-a2Ylim(2))*voxDim2v))*1e8)/1e8;
im3hStart = round((ax3Pos(1)/figPos(3)-((a3Xlim(1)-1)*voxDim3h))*1e8)/1e8;
im3hEnd = round(((ax3Pos(1)+ax3Pos(3))/figPos(3)+((handles.volume.imSz2-a3Xlim(2))*voxDim3h))*1e8)/1e8;
im3vStart = round((ax3Pos(2)/figPos(4)-((a3Ylim(1)-1)*voxDim3v))*1e8)/1e8;
im3vEnd = round(((ax3Pos(2)+ax3Pos(4))/figPos(4)+((handles.volume.imSz3-a3Ylim(2))*voxDim3v))*1e8)/1e8;

%axis arrow offset to exactly match crosshairs
s = settings;
try
    if s.matlab.desktop.HighDPISupport
        hAxOffset = 0.0005;
        vAxOffset = 0.0005;
    else
        hAxOffset = 0.0007;
        vAxOffset = 0.001;
    end
catch
    hAxOffset = 0.0007;
    vAxOffset = 0.001;
end


%PRECOMPUTE ARROWHEAD DYNAMIC LOCATIONS IN FIGURE NORMALIZED COORDINATES
if handles.volume.imSz1~=1
    handles.display.arrows.a1VarrowLocs = ...
        (im1hStart:(im1hEnd-im1hStart)/(handles.volume.imSz1-1):im1hEnd)-hAxOffset; %horizontal locations
    handles.display.arrows.a2VarrowLocs = ...
        (im2hStart:(im2hEnd-im2hStart)/(handles.volume.imSz1-1):im2hEnd)-hAxOffset; %horizontal locations
else
    handles.display.arrows.a1VarrowLocs = im1hStart; %horizontal location
    handles.display.arrows.a2VarrowLocs = im2hStart; %horizontal location
end
if handles.volume.imSz2~=1
    handles.display.arrows.a1HarrowLocs = ...
        (im1vStart:(im1vEnd-im1vStart)/(handles.volume.imSz2-1):im1vEnd)-vAxOffset; %vertical locations
    handles.display.arrows.a3VarrowLocs = ...
        (im3hStart:(im3hEnd-im3hStart)/(handles.volume.imSz2-1):im3hEnd)-hAxOffset; %horizontal locations
else
    handles.display.arrows.a1HarrowLocs = im1vStart; %vertical location
    handles.display.arrows.a3VarrowLocs = im3hStart; %horizontal location
end
if handles.volume.imSz3~=1
    handles.display.arrows.a2HarrowLocs = ...
        (im2vStart:(im2vEnd-im2vStart)/(handles.volume.imSz3-1):im2vEnd)-vAxOffset; %vertical locations
    handles.display.arrows.a3HarrowLocs = ...
        (im3vStart:(im3vEnd-im3vStart)/(handles.volume.imSz3-1):im3vEnd)-vAxOffset; %vertical locations
else
    handles.display.arrows.a2HarrowLocs = im2vStart; %vertical location
    handles.display.arrows.a3HarrowLocs = im3vStart; %vertical location
end

%DEFINE ARROWHEADS
arrowColor = [0.25 0.25 0.25];
arrowLength = 0.005*handles.figure.DPIscl;
headLength = 0.0062*figSz(1)*handles.figure.DPIscl;
headWidth = 0.0047*figSz(1)*handles.figure.DPIscl;

test = findall(0,'Tag','a1arrowH1');

if ~isempty(test) %ASSUME IF ONE ARROW FOUND, ALL EXIST
    %a1
    set(test,'X',[ax1hStart-arrowLength ax1hStart],'Y',...
        [handles.display.arrows.a1HarrowLocs(handles.volume.imSz2-handles.display.images.y+1) ...
        handles.display.arrows.a1HarrowLocs(handles.volume.imSz2-handles.display.images.y+1)],...
        'HeadLength',headLength,'HeadWidth',headWidth,'Color',arrowColor,'HitTest','off')
    handles.display.arrows.a1arrowH1 = test;
    test = findall(0,'Tag','a1arrowH2');
    set(test,'X',[ax1hEnd+arrowLength ax1hEnd],'Y',...
        [handles.display.arrows.a1HarrowLocs(handles.volume.imSz2-handles.display.images.y+1) ...
        handles.display.arrows.a1HarrowLocs(handles.volume.imSz2-handles.display.images.y+1)],...
        'HeadLength',headLength,'HeadWidth',headWidth,'Color',arrowColor,'HitTest','off')
    handles.display.arrows.a1arrowH2 = test;
    test = findall(0,'Tag','a1arrowV1');
    set(test,'X',[handles.display.arrows.a1VarrowLocs(handles.display.images.x) ...
        handles.display.arrows.a1VarrowLocs(handles.display.images.x)],...
        'Y',[ax1vStart-arrowLength ax1vStart],...
        'HeadLength',headLength,'HeadWidth',headWidth,'Color',arrowColor,'HitTest','off')
    handles.display.arrows.a1arrowV1 = test;
    test = findall(0,'Tag','a1arrowV2');
    set(test,'X',[handles.display.arrows.a1VarrowLocs(handles.display.images.x) ...
        handles.display.arrows.a1VarrowLocs(handles.display.images.x)],...
        'Y',[ax1vEnd+arrowLength ax1vEnd],...
        'HeadLength',headLength,'HeadWidth',headWidth,'Color',arrowColor,'HitTest','off')
    handles.display.arrows.a1arrowV2 = test;
    %a2
    test = findall(0,'Tag','a2arrowH1');
    set(test,'X',[ax2hStart-arrowLength ax2hStart],'Y',...
        [handles.display.arrows.a2HarrowLocs(handles.display.images.z) ...
        handles.display.arrows.a2HarrowLocs(handles.display.images.z)],...
        'HeadLength',headLength,'HeadWidth',headWidth,'Color',arrowColor,'HitTest','off')
    handles.display.arrows.a2arrowH1 = test;
    test = findall(0,'Tag','a2arrowH2');
    set(test,'X',[ax2hEnd+arrowLength ax2hEnd],'Y',...
        [handles.display.arrows.a2HarrowLocs(handles.display.images.z) ...
        handles.display.arrows.a2HarrowLocs(handles.display.images.z)],...
        'HeadLength',headLength,'HeadWidth',headWidth,'Color',arrowColor,'HitTest','off')
    handles.display.arrows.a2arrowH2 = test;
    test = findall(0,'Tag','a2arrowV1');
    set(test,'X',[handles.display.arrows.a2VarrowLocs(handles.display.images.x) ...
        handles.display.arrows.a2VarrowLocs(handles.display.images.x)],...
        'Y',[ax2vStart-arrowLength ax2vStart],...
        'HeadLength',headLength,'HeadWidth',headWidth,'Color',arrowColor,'HitTest','off')
    handles.display.arrows.a2arrowV1 = test;
    test = findall(0,'Tag','a2arrowV2');
    set(test,'X',[handles.display.arrows.a2VarrowLocs(handles.display.images.x) ...
        handles.display.arrows.a2VarrowLocs(handles.display.images.x)],...
        'Y',[ax2vEnd+arrowLength ax2vEnd],...
        'HeadLength',headLength,'HeadWidth',headWidth,'Color',arrowColor,'HitTest','off')
    handles.display.arrows.a2arrowV2 = test;
    %a3
    test = findall(0,'Tag','a3arrowH1');
    set(test,'X',[ax3hStart-arrowLength ax3hStart],'Y',...
        [handles.display.arrows.a3HarrowLocs(handles.display.images.z) ...
        handles.display.arrows.a3HarrowLocs(handles.display.images.z)],...
        'HeadLength',headLength,'HeadWidth',headWidth,'Color',arrowColor,'HitTest','off')
    handles.display.arrows.a3arrowH1 = test;
    test = findall(0,'Tag','a3arrowH2');
    set(test,'X',[ax3hEnd+arrowLength ax3hEnd],'Y',...
        [handles.display.arrows.a3HarrowLocs(handles.display.images.z) ...
        handles.display.arrows.a3HarrowLocs(handles.display.images.z)],...
        'HeadLength',headLength,'HeadWidth',headWidth,'Color',arrowColor,'HitTest','off')
    handles.display.arrows.a3arrowH2 = test;
    test = findall(0,'Tag','a3arrowV1');
    set(test,'X',[handles.display.arrows.a3VarrowLocs(handles.display.images.y) ...
        handles.display.arrows.a3VarrowLocs(handles.display.images.y)],...
        'Y',[ax3vStart-arrowLength ax3vStart],...
        'HeadLength',headLength,'HeadWidth',headWidth,'Color',arrowColor,'HitTest','off')
    handles.display.arrows.a3arrowV1 = test;
    test = findall(0,'Tag','a3arrowV2');
    set(test,'X',[handles.display.arrows.a3VarrowLocs(handles.display.images.y) ...
        handles.display.arrows.a3VarrowLocs(handles.display.images.y)],...
        'Y',[ax3vEnd+arrowLength ax3vEnd],...
        'HeadLength',headLength,'HeadWidth',headWidth,'Color',arrowColor,'HitTest','off')
    handles.display.arrows.a3arrowV2 = test;
else
    %a1
    handles.display.arrows.a1arrowH1 = ...
        annotation('arrow',[ax1hStart-arrowLength ax1hStart],...
        [handles.display.arrows.a1HarrowLocs(handles.volume.imSz2-handles.display.images.y+1) ...
        handles.display.arrows.a1HarrowLocs(handles.volume.imSz2-handles.display.images.y+1)],...
        'HeadLength',headLength,'HeadWidth',headWidth,'Color',arrowColor,'HitTest','off','Tag','a1arrowH1');
    handles.display.arrows.a1arrowH2 = ...
        annotation('arrow',[ax1hEnd+arrowLength ax1hEnd],...
        [handles.display.arrows.a1HarrowLocs(handles.volume.imSz2-handles.display.images.y+1) ...
        handles.display.arrows.a1HarrowLocs(handles.volume.imSz2-handles.display.images.y+1)],...
        'HeadLength',headLength,'HeadWidth',headWidth,'Color',arrowColor,'HitTest','off','Tag','a1arrowH2');
    handles.display.arrows.a1arrowV1 = ...
        annotation('arrow',[handles.display.arrows.a1VarrowLocs(handles.display.images.x) ...
        handles.display.arrows.a1VarrowLocs(handles.display.images.x)],...
        [ax1vStart-arrowLength ax1vStart],...
        'HeadLength',headLength,'HeadWidth',headWidth,'Color',arrowColor,'HitTest','off','Tag','a1arrowV1');
    handles.display.arrows.a1arrowV2 = ...
        annotation('arrow',[handles.display.arrows.a1VarrowLocs(handles.display.images.x) ...
        handles.display.arrows.a1VarrowLocs(handles.display.images.x)],...
        [ax1vEnd+arrowLength ax1vEnd],...
        'HeadLength',headLength,'HeadWidth',headWidth,'Color',arrowColor,'HitTest','off','Tag','a1arrowV2');
    %a2
    handles.display.arrows.a2arrowH1 = ...
        annotation('arrow',[ax2hStart-arrowLength ax2hStart],...
        [handles.display.arrows.a2HarrowLocs(handles.display.images.z) ...
        handles.display.arrows.a2HarrowLocs(handles.display.images.z)],...
        'HeadLength',headLength,'HeadWidth',headWidth,'Color',arrowColor,'HitTest','off','Tag','a2arrowH1');
    handles.display.arrows.a2arrowH2 = ...
        annotation('arrow',[ax2hEnd+arrowLength ax2hEnd],...
        [handles.display.arrows.a2HarrowLocs(handles.display.images.z) ...
        handles.display.arrows.a2HarrowLocs(handles.display.images.z)],...
        'HeadLength',headLength,'HeadWidth',headWidth,'Color',arrowColor,'HitTest','off','Tag','a2arrowH2');
    handles.display.arrows.a2arrowV1 = ...
        annotation('arrow',[handles.display.arrows.a2VarrowLocs(handles.display.images.x) ...
        handles.display.arrows.a2VarrowLocs(handles.display.images.x)],...
        [ax2vStart-arrowLength ax2vStart],...
        'HeadLength',headLength,'HeadWidth',headWidth,'Color',arrowColor,'HitTest','off','Tag','a2arrowV1');
    handles.display.arrows.a2arrowV2 = ...
        annotation('arrow',[handles.display.arrows.a2VarrowLocs(handles.display.images.x) ...
        handles.display.arrows.a2VarrowLocs(handles.display.images.x)],...
        [ax2vEnd+arrowLength ax2vEnd],...
        'HeadLength',headLength,'HeadWidth',headWidth,'Color',arrowColor,'HitTest','off','Tag','a2arrowV2');
    %a3
    handles.display.arrows.a3arrowH1 = ...
        annotation('arrow',[ax3hStart-arrowLength ax3hStart],...
        [handles.display.arrows.a3HarrowLocs(handles.display.images.z) ...
        handles.display.arrows.a3HarrowLocs(handles.display.images.z)],...
        'HeadLength',headLength,'HeadWidth',headWidth,'Color',arrowColor,'HitTest','off','Tag','a3arrowH1');
    handles.display.arrows.a3arrowH2 = ...
        annotation('arrow',[ax3hEnd+arrowLength ax3hEnd],...
        [handles.display.arrows.a3HarrowLocs(handles.display.images.z) ...
        handles.display.arrows.a3HarrowLocs(handles.display.images.z)],...
        'HeadLength',headLength,'HeadWidth',headWidth,'Color',arrowColor,'HitTest','off','Tag','a3arrowH2');
    handles.display.arrows.a3arrowV1 = ...
        annotation('arrow',[handles.display.arrows.a3VarrowLocs(handles.display.images.y) ...
        handles.display.arrows.a3VarrowLocs(handles.display.images.y)],...
        [ax3vStart-arrowLength ax3vStart],...
        'HeadLength',headLength,'HeadWidth',headWidth,'Color',arrowColor,'HitTest','off','Tag','a3arrowV1');
    handles.display.arrows.a3arrowV2 = ...
        annotation('arrow',[handles.display.arrows.a3VarrowLocs(handles.display.images.y) ...
        handles.display.arrows.a3VarrowLocs(handles.display.images.y)],...
        [ax3vEnd+arrowLength ax3vEnd],...
        'HeadLength',headLength,'HeadWidth',headWidth,'Color',arrowColor,'HitTest','off','Tag','a3arrowV2');
end

end

function maxDispSliderCallback(hObj,~)

global handles

val = getSliderValue(hObj,handles.display.images.sclFctr);

while val~=handles.display.images.cLim(2)
    
    set(hObj,'Value',val)
    handles.display.images.cLim(2) = val;
    
    if handles.display.images.cLim(1)>=val
        stepDist = abs(handles.display.sliders.s4steps-val);
        idx = find(stepDist==min(stepDist));
        set(handles.display.sliders.s5,'Value',handles.display.sliders.s5steps(idx))
        handles.display.images.cLim(1) = handles.display.sliders.s5steps(idx);
    end
    
    s4tString = formatValueString(handles.display.images.cLim(2)/handles.display.images.sclFctr);
    s5tString = formatValueString(handles.display.images.cLim(1)/handles.display.images.sclFctr);
    
    set(handles.display.sliders.s4t,'String',s4tString)
    set(handles.display.sliders.s5t,'String',s5tString)
    
    if strcmpi(handles.volume.modality,'CT')
        test = get(findobj('Tag','CMPanel'),'Children');
        if ~isempty(test), set(test,'Value',0); end
        if all(handles.display.images.cLim==getCTwindow('Body'))
            set(test(strcmpi(get(test,'String'),'Body')),'Value',1)
        elseif all(handles.display.images.cLim==getCTwindow('Lung'))
            set(test(strcmpi(get(test,'String'),'Lung')),'Value',1)
        elseif all(handles.display.images.cLim==getCTwindow('Bone'))
            set(test(strcmpi(get(test,'String'),'Bone')),'Value',1)
        end
    end
    
    drawnow
    
    setDisplayWindow()
    setDisplaySliders()
    
    drawnow
    
    val = getSliderValue(hObj,handles.display.images.sclFctr);
    
end

drawnow

set(handles.display.sliders.s4,'Value',handles.display.images.cLim(2))
set(handles.display.sliders.s5,'Value',handles.display.images.cLim(1))

end

function minDispSliderCallback(hObj,~)

global handles

val = getSliderValue(hObj,handles.display.images.sclFctr);

while val~=handles.display.images.cLim(1)
    
    set(hObj,'Value',val)
    handles.display.images.cLim(1) = val;
    
    if val>=handles.display.images.cLim(2)
        stepDist = abs(handles.display.sliders.s5steps-val);
        idx = find(stepDist==min(stepDist));
        set(handles.display.sliders.s4,'Value',handles.display.sliders.s4steps(idx))
        handles.display.images.cLim(2) = handles.display.sliders.s4steps(idx);
    end
    
    s4tString = formatValueString(handles.display.images.cLim(2)/handles.display.images.sclFctr);
    s5tString = formatValueString(handles.display.images.cLim(1)/handles.display.images.sclFctr);
    
    set(handles.display.sliders.s4t,'String',s4tString)
    set(handles.display.sliders.s5t,'String',s5tString)
    
    if strcmpi(handles.volume.modality,'CT')
        test = get(findobj('Tag','CMPanel'),'Children');
        if ~isempty(test), set(test,'Value',0); end
        if all(handles.display.images.cLim==getCTwindow('Body'))
            set(test(strcmpi(get(test,'String'),'Body')),'Value',1)
        elseif all(handles.display.images.cLim==getCTwindow('Lung'))
            set(test(strcmpi(get(test,'String'),'Lung')),'Value',1)
        elseif all(handles.display.images.cLim==getCTwindow('Bone'))
            set(test(strcmpi(get(test,'String'),'Bone')),'Value',1)
        end
    end
    
    drawnow
    
    setDisplayWindow()
    setDisplaySliders()
    
    drawnow
    
    val = getSliderValue(hObj,handles.display.images.sclFctr);
    
end

drawnow

set(handles.display.sliders.s4,'Value',handles.display.images.cLim(2))
set(handles.display.sliders.s5,'Value',handles.display.images.cLim(1))

end

function maxDispSliderCallback_1(hObj,~)

global handles

val = getSliderValue(hObj,handles.display.images.sclFctr_1);

while val~=handles.display.images.cLim_1(2)
    
    set(hObj,'Value',val)
    handles.display.images.cLim_1(2) = val;
    
    if handles.display.images.cLim_1(1)>=val
        stepDist = abs(handles.display.sliders.s4steps_1-val);
        idx = find(stepDist==min(stepDist));
        set(handles.display.sliders.s5_1,'Value',handles.display.sliders.s5steps_1(idx))
        handles.display.images.cLim_1(1) = handles.display.sliders.s5steps_1(idx);
    end
    
    s4tString = formatValueString(handles.display.images.cLim_1(2)/handles.display.images.sclFctr_1);
    s5tString = formatValueString(handles.display.images.cLim_1(1)/handles.display.images.sclFctr_1);
    
    set(handles.display.sliders.s4t_1,'String',s4tString)
    set(handles.display.sliders.s5t_1,'String',s5tString)
    
    if strcmpi(handles.volume.modality_1,'CT')
        test = get(findobj('Tag','CMPanel_1'),'Children');
        if ~isempty(test), set(test,'Value',0); end
        if all(handles.display.images.cLim_1==getCTwindow('Body'))
            set(test(strcmpi(get(test,'String'),'Body')),'Value',1)
        elseif all(handles.display.images.cLim_1==getCTwindow('Lung'))
            set(test(strcmpi(get(test,'String'),'Lung')),'Value',1)
        elseif all(handles.display.images.cLim_1==getCTwindow('Bone'))
            set(test(strcmpi(get(test,'String'),'Bone')),'Value',1)
        end
    end
    
    drawnow
    
    setFusedDisplayWindow()
    setDisplaySliders_1()
    
    drawnow
    
    val = getSliderValue(hObj,handles.display.images.sclFctr_1);
    
end

set(handles.display.sliders.s4_1,'Value',handles.display.images.cLim_1(2))
set(handles.display.sliders.s5_1,'Value',handles.display.images.cLim_1(1))

drawnow

end

function minDispSliderCallback_1(hObj,~)

global handles

val = getSliderValue(hObj,handles.display.images.sclFctr_1);

while val~=handles.display.images.cLim_1(1)
    
    set(hObj,'Value',val)
    handles.display.images.cLim_1(1) = val;
    
    if val>=handles.display.images.cLim_1(2)
        stepDist = abs(handles.display.sliders.s5steps_1-val);
        idx = find(stepDist==min(stepDist));
        set(handles.display.sliders.s4_1,'Value',handles.display.sliders.s4steps_1(idx))
        handles.display.images.cLim_1(2) = handles.display.sliders.s4steps_1(idx);
    end
    
    s4tString = formatValueString(handles.display.images.cLim_1(2)/handles.display.images.sclFctr_1);
    s5tString = formatValueString(handles.display.images.cLim_1(1)/handles.display.images.sclFctr_1);
    
    set(handles.display.sliders.s4t_1,'String',s4tString)
    set(handles.display.sliders.s5t_1,'String',s5tString)
    
    if strcmpi(handles.volume.modality_1,'CT')
        test = get(findobj('Tag','CMPanel_1'),'Children');
        if ~isempty(test), set(test,'Value',0); end
        if all(handles.display.images.cLim_1==getCTwindow('Body'))
            set(test(strcmpi(get(test,'String'),'Body')),'Value',1)
        elseif all(handles.display.images.cLim_1==getCTwindow('Lung'))
            set(test(strcmpi(get(test,'String'),'Lung')),'Value',1)
        elseif all(handles.display.images.cLim_1==getCTwindow('Bone'))
            set(test(strcmpi(get(test,'String'),'Bone')),'Value',1)
        end
    end
    
    drawnow
    
    setFusedDisplayWindow()
    setDisplaySliders_1()
    
    drawnow
    
    val = getSliderValue(hObj,handles.display.images.sclFctr_1);
    
end

set(handles.display.sliders.s4_1,'Value',handles.display.images.cLim_1(2))
set(handles.display.sliders.s5_1,'Value',handles.display.images.cLim_1(1))

drawnow

end

function maxDispSliderCallback_2(hObj,~)

global handles

val = getSliderValue(hObj,handles.display.images.sclFctr_2);

while val~=handles.display.images.cLim_2(2)
    
    set(hObj,'Value',val)
    handles.display.images.cLim_2(2) = val;
    
    if handles.display.images.cLim_2(1)>=val
        stepDist = abs(handles.display.sliders.s4steps_2-val);
        idx = find(stepDist==min(stepDist));
        set(handles.display.sliders.s5_2,'Value',handles.display.sliders.s5steps_2(idx))
        handles.display.images.cLim_2(1) = handles.display.sliders.s5steps_2(idx);
    end
    
    s4tString = formatValueString(handles.display.images.cLim_2(2)/handles.display.images.sclFctr_2);
    s5tString = formatValueString(handles.display.images.cLim_2(1)/handles.display.images.sclFctr_2);
    
    set(handles.display.sliders.s4t_2,'String',s4tString)
    set(handles.display.sliders.s5t_2,'String',s5tString)
    
    if strcmpi(handles.volume.modality_2,'CT')
        test = get(findobj('Tag','CMPanel_2'),'Children');
        if ~isempty(test), set(test,'Value',0); end
        if all(handles.display.images.cLim_2==getCTwindow('Body'))
            set(test(strcmpi(get(test,'String'),'Body')),'Value',1)
        elseif all(handles.display.images.cLim_2==getCTwindow('Lung'))
            set(test(strcmpi(get(test,'String'),'Lung')),'Value',1)
        elseif all(handles.display.images.cLim_2==getCTwindow('Bone'))
            set(test(strcmpi(get(test,'String'),'Bone')),'Value',1)
        end
    end
    
    drawnow
    
    setFusedDisplayWindow()
    setDisplaySliders_2()
    
    drawnow
    
    val = getSliderValue(hObj,handles.display.images.sclFctr_2);
    
end

set(handles.display.sliders.s4_2,'Value',handles.display.images.cLim_2(2))
set(handles.display.sliders.s5_2,'Value',handles.display.images.cLim_2(1))

drawnow

end

function minDispSliderCallback_2(hObj,~)

global handles

val = getSliderValue(hObj,handles.display.images.sclFctr_2);

while val~=handles.display.images.cLim_2(1)
    
    set(hObj,'Value',val)
    handles.display.images.cLim_2(1) = val;
    
    if val>=handles.display.images.cLim_2(2)
        stepDist = abs(handles.display.sliders.s5steps_2-val);
        idx = find(stepDist==min(stepDist));
        set(handles.display.sliders.s4_2,'Value',handles.display.sliders.s4steps_2(idx))
        handles.display.images.cLim_2(2) = handles.display.sliders.s4steps_2(idx);
    end
    
    s4tString = formatValueString(handles.display.images.cLim_2(2)/handles.display.images.sclFctr_2);
    s5tString = formatValueString(handles.display.images.cLim_2(1)/handles.display.images.sclFctr_2);
    
    set(handles.display.sliders.s4t_2,'String',s4tString)
    set(handles.display.sliders.s5t_2,'String',s5tString)
    
    if strcmpi(handles.volume.modality_2,'CT')
        test = get(findobj('Tag','CMPanel_2'),'Children');
        if ~isempty(test), set(test,'Value',0); end
        if all(handles.display.images.cLim_2==getCTwindow('Body'))
            set(test(strcmpi(get(test,'String'),'Body')),'Value',1)
        elseif all(handles.display.images.cLim_2==getCTwindow('Lung'))
            set(test(strcmpi(get(test,'String'),'Lung')),'Value',1)
        elseif all(handles.display.images.cLim_2==getCTwindow('Bone'))
            set(test(strcmpi(get(test,'String'),'Bone')),'Value',1)
        end
    end
    
    drawnow
    
    setFusedDisplayWindow()
    setDisplaySliders_2()
    
    drawnow
    
    val = getSliderValue(hObj,handles.display.images.sclFctr_2);
    
end

set(handles.display.sliders.s4_2,'Value',handles.display.images.cLim_2(2))
set(handles.display.sliders.s5_2,'Value',handles.display.images.cLim_2(1))

drawnow

end

function value = getSliderValue(hObj,sclFctr)

int = 100;

value = single(get(hObj,'Value'));
if (get(hObj,'Max')-get(hObj,'Min'))/sclFctr>int
    value = round(value);
end
%account for casting errors
if value>get(hObj,'Max')
    value = get(hObj,'Max');
end
if value<get(hObj,'Min')
    value = get(hObj,'Min');
end

end

function maxDispTextCallback(hObj,~)

global handles

int = 100;

val = str2double(get(hObj,'String'));

%assure valid input
if (isempty(val) || ~isfinite(val) || ~isreal(val))
    val = get(handles.display.sliders.s4,'Value');
end
if val<get(handles.display.sliders.s4,'Min')/handles.display.images.sclFctr
    val = get(handles.display.sliders.s4,'Min')/handles.display.images.sclFctr;
end
if val>get(handles.display.sliders.s4,'Max')/handles.display.images.sclFctr
    val = get(handles.display.sliders.s4,'Max')/handles.display.images.sclFctr;
end

if (get(handles.display.sliders.s4,'Max')-get(handles.display.sliders.s4,'Min'))/handles.display.images.sclFctr>int
    val = round(val);
end

str = formatValueString(val);
set(hObj,'String',str)

set(handles.display.sliders.s4,'Value',val*handles.display.images.sclFctr)
maxDispSliderCallback(handles.display.sliders.s4)

end

function minDispTextCallback(hObj,~)

global handles

int = 100;

val = str2double(get(hObj,'String'));

%assure valid input
if (isempty(val) || ~isfinite(val) || ~isreal(val))
    val = get(handles.display.sliders.s5,'Value');
end
if val<get(handles.display.sliders.s5,'Min')/handles.display.images.sclFctr
    val = get(handles.display.sliders.s5,'Min')/handles.display.images.sclFctr;
end
if val>get(handles.display.sliders.s5,'Max')/handles.display.images.sclFctr
    val = get(handles.display.sliders.s5,'Max')/handles.display.images.sclFctr;
end

if (get(handles.display.sliders.s4,'Max')-get(handles.display.sliders.s4,'Min'))/handles.display.images.sclFctr>int
    val = round(val);
end

str = formatValueString(val);
set(hObj,'String',str)

set(handles.display.sliders.s5,'Value',val*handles.display.images.sclFctr)
minDispSliderCallback(handles.display.sliders.s5)

end

function maxDispTextCallback_1(hObj,~)

global handles

int = 100;

val = str2double(get(hObj,'String'));

%assure valid input
if (isempty(val) || ~isfinite(val) || ~isreal(val))
    val = get(handles.display.sliders.s4_1,'Value');
end
if val<get(handles.display.sliders.s4_1,'Min')/handles.display.images.sclFctr_1
    val = get(handles.display.sliders.s4_1,'Min')/handles.display.images.sclFctr_1;
end
if val>get(handles.display.sliders.s4_1,'Max')/handles.display.images.sclFctr_1
    val = get(handles.display.sliders.s4_1,'Max')/handles.display.images.sclFctr_1;
end

if (get(handles.display.sliders.s4_1,'Max')-get(handles.display.sliders.s4_1,'Min'))/handles.display.images.sclFctr_1>int
    val = round(val);
end

str = formatValueString(val);
set(hObj,'String',str)

set(handles.display.sliders.s4_1,'Value',val*handles.display.images.sclFctr_1)
maxDispSliderCallback_1(handles.display.sliders.s4_1)

end

function minDispTextCallback_1(hObj,~)

global handles

int = 100;

val = str2double(get(hObj,'String'));

%assure valid input
if (isempty(val) || ~isfinite(val) || ~isreal(val))
    val = get(handles.display.sliders.s5_1,'Value');
end
if val<get(handles.display.sliders.s5_1,'Min')/handles.display.images.sclFctr_1
    val = get(handles.display.sliders.s5_1,'Min')/handles.display.images.sclFctr_1;
end
if val>get(handles.display.sliders.s5_1,'Max')/handles.display.images.sclFctr_1
    val = get(handles.display.sliders.s5_1,'Max')/handles.display.images.sclFctr_1;
end

if (get(handles.display.sliders.s5_1,'Max')-get(handles.display.sliders.s5_1,'Min'))/handles.display.images.sclFctr_1>int
    val = round(val);
end

str = formatValueString(val);
set(hObj,'String',str)

set(handles.display.sliders.s5_1,'Value',val*handles.display.images.sclFctr_1)
minDispSliderCallback_1(handles.display.sliders.s5_1)

end

function maxDispTextCallback_2(hObj,~)

global handles

int = 100;

val = str2double(get(hObj,'String'));

%assure valid input
if (isempty(val) || ~isfinite(val) || ~isreal(val))
    val = get(handles.display.sliders.s4_2,'Value');
end
if val<get(handles.display.sliders.s4_2,'Min')/handles.display.images.sclFctr_2
    val = get(handles.display.sliders.s4_2,'Min')/handles.display.images.sclFctr_2;
end
if val>get(handles.display.sliders.s4_2,'Max')/handles.display.images.sclFctr_2
    val = get(handles.display.sliders.s4_2,'Max')/handles.display.images.sclFctr_2;
end

if (get(handles.display.sliders.s4_2,'Max')-get(handles.display.sliders.s4_2,'Min'))/handles.display.images.sclFctr_2>int
    val = round(val);
end

str = formatValueString(val);
set(hObj,'String',str)

set(handles.display.sliders.s4_2,'Value',val*handles.display.images.sclFctr_2)
maxDispSliderCallback_2(handles.display.sliders.s4_2)

end

function minDispTextCallback_2(hObj,~)

global handles

int = 100;

val = str2double(get(hObj,'String'));

%assure valid input
if (isempty(val) || ~isfinite(val) || ~isreal(val))
    val = get(handles.display.sliders.s5_2,'Value');
end
if val<get(handles.display.sliders.s5_2,'Min')/handles.display.images.sclFctr_2
    val = get(handles.display.sliders.s5_2,'Min')/handles.display.images.sclFctr_2;
end
if val>get(handles.display.sliders.s5_2,'Max')/handles.display.images.sclFctr_2
    val = get(handles.display.sliders.s5_2,'Max')/handles.display.images.sclFctr_2;
end

if (get(handles.display.sliders.s5_2,'Max')-get(handles.display.sliders.s5_2,'Min'))/handles.display.images.sclFctr_2>int
    val = round(val);
end

str = formatValueString(val);
set(hObj,'String',str)

set(handles.display.sliders.s5_2,'Value',val*handles.display.images.sclFctr_2)
minDispSliderCallback_2(handles.display.sliders.s5_2)

end

function maxDispButtonCallback(hObj,~)

global handles

%find max coordinates
idx = find(handles.volume.vol==handles.display.images.volMax);
idx = idx(1);
[x,y,z] = ind2sub([handles.volume.imSz1,handles.volume.imSz2,handles.volume.imSz3],idx);

updateCurrentVoxel(x,y,z)

end

function minDispButtonCallback(hObj,~)

global handles

%find min coordinates
idx = find(handles.volume.vol==handles.display.images.volMin);
idx = idx(1);
[x,y,z] = ind2sub([handles.volume.imSz1,handles.volume.imSz2,handles.volume.imSz3],idx);

updateCurrentVoxel(x,y,z)

end

function maxDispButtonCallback_1(hObj,~)

global handles

%find max coordinates
idx = find(handles.volume.vol_1==handles.display.images.volMax_1);
idx = idx(1);
[x,y,z] = ind2sub([handles.volume.imSz1,handles.volume.imSz2,handles.volume.imSz3],idx);

updateCurrentVoxel(x,y,z)

end

function minDispButtonCallback_1(hObj,~)

global handles

%find min coordinates
idx = find(handles.volume.vol_1==handles.display.images.volMin_1);
idx = idx(1);
[x,y,z] = ind2sub([handles.volume.imSz1,handles.volume.imSz2,handles.volume.imSz3],idx);

updateCurrentVoxel(x,y,z)

end

function maxDispButtonCallback_2(hObj,~)

global handles

%find max coordinates
idx = find(handles.volume.vol_2==handles.display.images.volMax_2);
idx = idx(1);
[x,y,z] = ind2sub([handles.volume.imSz1,handles.volume.imSz2,handles.volume.imSz3],idx);

updateCurrentVoxel(x,y,z)

end

function minDispButtonCallback_2(hObj,~)

global handles

%find min coordinates
idx = find(handles.volume.vol_2==handles.display.images.volMin_2);
idx = idx(1);
[x,y,z] = ind2sub([handles.volume.imSz1,handles.volume.imSz2,handles.volume.imSz3],idx);

updateCurrentVoxel(x,y,z)

end

function updateCurrentVoxel(x,y,z)

global handles

handles.display.images.x = x;
handles.display.images.y = y;
handles.display.images.z = z;

if ~handles.display.images.fusedImages
    %images
    showImage1(handles.volume.vol,handles.display.images.h1,handles.display.images.z,handles.display.images.fltr)
    showImage2(handles.volume.vol,handles.display.images.h2,handles.display.images.y,handles.display.images.fltr)
    showImage3(handles.volume.vol,handles.display.images.h3,handles.display.images.x,handles.display.images.fltr)
    %sliders
    set(handles.display.sliders.s1,'Value',handles.display.images.z)
    set(handles.display.sliders.s2,'Value',handles.display.images.y)
    set(handles.display.sliders.s3,'Value',handles.display.images.x)
else
    %images
    showFusedImage1(handles.axes.a1,handles.volume.vol_1,handles.volume.vol_2,...
        handles.display.images.h1_1,handles.display.images.h1_2,...
        handles.display.images.z,handles.display.images.cLim_1,handles.display.images.cLim_2,...
        handles.display.images.cm_1,handles.display.images.cm_2,handles.display.images.fltr)
    showFusedImage2(handles.axes.a2,handles.volume.vol_1,handles.volume.vol_2,...
        handles.display.images.h2_1,handles.display.images.h2_2,...
        handles.display.images.y,handles.display.images.cLim_1,handles.display.images.cLim_2,...
        handles.display.images.cm_1,handles.display.images.cm_2,handles.display.images.fltr)
    showFusedImage3(handles.axes.a3,handles.volume.vol_1,handles.volume.vol_2,...
        handles.display.images.h3_1,handles.display.images.h3_2,...
        handles.display.images.x,handles.display.images.cLim_1,handles.display.images.cLim_2,...
        handles.display.images.cm_1,handles.display.images.cm_2,handles.display.images.fltr)
    %sliders
    set(handles.display.sliders.s1,'Value',handles.display.images.z)
    set(handles.display.sliders.s2,'Value',handles.display.images.y)
    set(handles.display.sliders.s3,'Value',handles.display.images.x)
end

if handles.display.images.drawROI
    showMask()
end

updateZoomWindow

end

function smSliderCallback(hObj,~)

global handles

FWHM = round(get(hObj,'Value'));
s1str = get(handles.display.sliders.s1t,'String');
s2str = get(handles.display.sliders.s2t,'String');
s3str = get(handles.display.sliders.s3t,'String');
s4str = get(handles.display.sliders.s4t,'String');
s5str = get(handles.display.sliders.s5t,'String');

set(hObj,'Value',FWHM)
set(handles.display.sliders.s6t,'String',FWHM)

if eval(['isfield(handles.volume,''vol' num2str(FWHM) 'mm'')'])
    eval(['handles.volume.vol = handles.volume.vol' num2str(FWHM) 'mm;'])
else
    set(handles.figure.f1,'Pointer','watch')
    makeInactive()
    
    %SMOOTHING TEXT
    set(findobj('Tag','s6Text2'),'Visible','on')
    %VOLUME STAMP TEXT
    set(findobj('Tag','StampText1'),'Visible','on')
    set(findobj('Tag','StampText2'),'Visible','on')
    drawnow
    
    %ACCOUNT FOR ASYMMETRIC VOXELS
    siz = FWHM./[handles.volume.voxDimX handles.volume.voxDimY handles.volume.voxDimZ];
    %SMOOTH VOLUME
    eval(['handles.volume.vol' num2str(FWHM) 'mm = GaussianSmooth(handles.volume.vol0mm, siz);'])
    eval(['handles.volume.vol = handles.volume.vol' num2str(FWHM) 'mm;'])
    
    set(handles.figure.f1,'Pointer','arrow')
    makeActive()
end

set(handles.display.sliders.s1t,'String',s1str)
set(handles.display.sliders.s2t,'String',s2str)
set(handles.display.sliders.s3t,'String',s3str)
set(handles.display.sliders.s4t,'String',s4str)
set(handles.display.sliders.s5t,'String',s5str)
set(handles.display.sliders.s6t,'String',FWHM)

handles.display.images.h1 = setVolume(findobj('Tag','a1im'),handles.axes.a1);
handles.display.images.h2 = setVolume(findobj('Tag','a2im'),handles.axes.a2);
handles.display.images.h3 = setVolume(findobj('Tag','a3im'),handles.axes.a3);

if ~handles.display.images.fusedImages
    figMotionFcn()
else
    figMotionFcnFused()
end

end

function smSliderText(hObj,~)

global handles

set(handles.display.sliders.s6t,'String',num2str(round(get(hObj,'Value'))))

end

function smTextCallback(hObj,~)

global handles

val = str2double(get(hObj,'String'));

%assure valid input
if (isempty(val) || ~isfinite(val) || ~isreal(val)) || ...
        val<get(handles.display.sliders.s6,'Min') || val>get(handles.display.sliders.s6,'Max')
    val = get(handles.display.sliders.s6,'Value');
    set(hObj,'String',num2str(val))
end

set(handles.display.sliders.s6,'Value',val)
smSliderCallback(handles.display.sliders.s6)

end

function fwSliderCallback(hObj,~)

global handles

alph = get(hObj,'Value');

alpha(handles.display.images.h1_1,1-alph);
alpha(handles.display.images.h2_1,1-alph);
alpha(handles.display.images.h3_1,1-alph);

end

function a1SliderListener(hObj,~)

global handles

holdSlider()

handles.display.images.z = round(get(hObj,'Value'));
if ~handles.display.images.fusedImages
    showImage1(handles.volume.vol,handles.display.images.h1,handles.display.images.z,handles.display.images.fltr)
else
    showFusedImage1(handles.axes.a1,handles.volume.vol_1,handles.volume.vol_2,...
        handles.display.images.h1_1,handles.display.images.h1_2,...
        handles.display.images.z,handles.display.images.cLim_1,handles.display.images.cLim_2,...
        handles.display.images.cm_1,handles.display.images.cm_2,handles.display.images.fltr)
end
moveArrows()
%SLIDER TEXT
axesSliderText(handles.display.sliders.s1t);

if handles.display.images.drawROI
    showMask()
end

end

function a2SliderListener(hObj,~)

global handles

holdSlider()

handles.display.images.y = round(get(hObj,'Value'));
if ~handles.display.images.fusedImages
    showImage2(handles.volume.vol,handles.display.images.h2,handles.display.images.y,handles.display.images.fltr)
else
    showFusedImage2(handles.axes.a2,handles.volume.vol_1,handles.volume.vol_2,...
        handles.display.images.h2_1,handles.display.images.h2_2,...
        handles.display.images.y,handles.display.images.cLim_1,handles.display.images.cLim_2,...
        handles.display.images.cm_1,handles.display.images.cm_2,handles.display.images.fltr)
end
moveArrows()
%SLIDER TEXT
axesSliderText(handles.display.sliders.s2t);

if handles.display.images.drawROI
    showMask()
end

end

function a3SliderListener(hObj,~)

global handles

holdSlider()

handles.display.images.x = round(get(hObj,'Value'));
if ~handles.display.images.fusedImages
    showImage3(handles.volume.vol,handles.display.images.h3,handles.display.images.x,handles.display.images.fltr)
else
    showFusedImage3(handles.axes.a3,handles.volume.vol_1,handles.volume.vol_2,...
        handles.display.images.h3_1,handles.display.images.h3_2,...
        handles.display.images.x,handles.display.images.cLim_1,handles.display.images.cLim_2,...
        handles.display.images.cm_1,handles.display.images.cm_2,handles.display.images.fltr)
end
moveArrows()
%SLIDER TEXT
axesSliderText(handles.display.sliders.s3t);

if handles.display.images.drawROI
    showMask()
end

end

function a1TextCallback(hObj,~)

global handles

val = str2double(get(hObj,'String'));

%assure valid input
if (isempty(val) || ~isfinite(val) || ~isreal(val)) || ...
        val<get(handles.display.sliders.s1,'Min') || val>get(handles.display.sliders.s1,'Max')
    val = get(handles.display.sliders.s1,'Value');
    set(hObj,'String',num2str(val))
end

set(handles.display.sliders.s1,'Value',val)
a1SliderListener(handles.display.sliders.s1)

releaseSlider()

end

function a2TextCallback(hObj,~)

global handles

val = str2double(get(hObj,'String'));

%assure valid input
if (isempty(val) || ~isfinite(val) || ~isreal(val)) || ...
        val<get(handles.display.sliders.s2,'Min') || val>get(handles.display.sliders.s2,'Max')
    val = get(handles.display.sliders.s2,'Value');
    set(hObj,'String',num2str(val))
end

set(handles.display.sliders.s2,'Value',val)
a2SliderListener(handles.display.sliders.s2)

releaseSlider()

end

function a3TextCallback(hObj,~)

global handles

val = str2double(get(hObj,'String'));

%assure valid input
if (isempty(val) || ~isfinite(val) || ~isreal(val)) || ...
        val<get(handles.display.sliders.s3,'Min') || val>get(handles.display.sliders.s3,'Max')
    val = get(handles.display.sliders.s3,'Value');
    set(hObj,'String',num2str(val))
end

set(handles.display.sliders.s3,'Value',val)
a3SliderListener(handles.display.sliders.s3)

releaseSlider()

end

function holdSlider

global handles

handles.display.sliders.sliderGrab = true;

end

function releaseSlider(varargin)

global handles

handles.display.sliders.sliderGrab = false;
if ~handles.display.images.fusedImages
    figMotionFcn()
else
    figMotionFcnFused()
end

end

function figMotionFcn(varargin)

global handles

%ASSUME EVERY UI ELEMENT DEFINED IN NORMALIZED VALUES
units = get(handles.figure.f1,'Units');
set(handles.figure.f1,'Units','normalized');

pt = get(handles.figure.f1,'CurrentPoint');
x = pt(1,1); y = pt(1,2);

set(handles.figure.f1,'Units',units);

if x>=handles.axes.a1Pos(1) && ...
        x<(handles.axes.a1Pos(1)+handles.axes.a1Pos(3)) && ...
        y>=handles.axes.a1Pos(2) && ...
        y<(handles.axes.a1Pos(2)+handles.axes.a1Pos(4))
    if ~handles.display.sliders.sliderGrab
        set(handles.figure.f1,'Pointer','crosshair')
        pt = get(handles.axes.a1,'CurrentPoint');
        pt_x = round(pt(1,1)); pt_y = round(pt(1,2));
        [x,y,z] = getCoordinates(handles.axes.a1,pt_x,pt_y);
        coordinatesValue(x,y,z)
        voxelValue(x,y,z)
        if handles.display.text.showLocs
            locValue(x,y,z)
        end
        turnOnText()
    else
        set(handles.figure.f1,'Pointer','hand')
        turnOffText()
        turnOffCrosshair()
    end
    moveCrosshair(handles.axes.a1)
    handles.display.cursor.overAxes = true;
elseif x>=handles.axes.a2Pos(1) && ...
        x<(handles.axes.a2Pos(1)+handles.axes.a2Pos(3)) && ...
        y>=handles.axes.a2Pos(2) && ...
        y<(handles.axes.a2Pos(2)+handles.axes.a2Pos(4))
    if ~handles.display.sliders.sliderGrab
        set(handles.figure.f1,'Pointer','crosshair')
        pt = get(handles.axes.a2,'CurrentPoint');
        pt_x = round(pt(1,1)); pt_y = handles.volume.imSz3-round(pt(1,2))+1;
        [x,y,z] = getCoordinates(handles.axes.a2,pt_x,pt_y);
        coordinatesValue(x,y,z)
        voxelValue(x,y,z)
        if handles.display.text.showLocs
            locValue(x,y,z)
        end
        turnOnText()
    else
        set(handles.figure.f1,'Pointer','hand')
        turnOffText()
        turnOffCrosshair()
    end
    moveCrosshair(handles.axes.a2)
    handles.display.cursor.overAxes = true;
elseif x>=handles.axes.a3Pos(1) && ...
        x<(handles.axes.a3Pos(1)+handles.axes.a3Pos(3)) && ...
        y>=handles.axes.a3Pos(2) && ...
        y<(handles.axes.a3Pos(2)+handles.axes.a3Pos(4))
    if ~handles.display.sliders.sliderGrab
        set(handles.figure.f1,'Pointer','crosshair')
        pt = get(handles.axes.a3,'CurrentPoint');
        pt_x = round(pt(1,1)); pt_y = handles.volume.imSz3-round(pt(1,2))+1;
        [x,y,z] = getCoordinates(handles.axes.a3,pt_x,pt_y);
        coordinatesValue(x,y,z)
        voxelValue(x,y,z)
        if handles.display.text.showLocs
            locValue(x,y,z)
        end
        turnOnText()
    else
        set(handles.figure.f1,'Pointer','hand')
        turnOffText()
        turnOffCrosshair()
    end
    moveCrosshair(handles.axes.a3)
    handles.display.cursor.overAxes = true;
elseif x>=handles.display.sliders.s1Pos(1) && ...
        x<handles.display.sliders.s1Pos(1)+handles.display.sliders.s1Pos(3) ...
        && y>=handles.display.sliders.s1Pos(2) && ...
        y<handles.display.sliders.s1Pos(2)+handles.display.sliders.s1Pos(4)
    if ~handles.display.sliders.sliderGrab
        if strcmpi(get(handles.display.sliders.s1,'Enable'),'on')
            set(handles.figure.f1,'Pointer','hand')
        else
            set(handles.figure.f1,'Pointer','arrow')
        end
    else
        set(handles.figure.f1,'Pointer','hand')
    end
    turnOffText()
    turnOffCrosshair()
    handles.display.cursor.overAxes = false;
elseif x>=handles.display.sliders.s2Pos(1) && ...
        x<handles.display.sliders.s2Pos(1)+handles.display.sliders.s2Pos(3) && ...
        y>=handles.display.sliders.s2Pos(2) && ...
        y<handles.display.sliders.s2Pos(2)+handles.display.sliders.s2Pos(4)
    if ~handles.display.sliders.sliderGrab
        if strcmpi(get(handles.display.sliders.s2,'Enable'),'on')
            set(handles.figure.f1,'Pointer','hand')
        else
            set(handles.figure.f1,'Pointer','arrow')
        end
    else
        set(handles.figure.f1,'Pointer','hand')
    end
    turnOffText()
    turnOffCrosshair()
    handles.display.cursor.overAxes = false;
elseif x>=handles.display.sliders.s3Pos(1) && ...
        x<handles.display.sliders.s3Pos(1)+handles.display.sliders.s3Pos(3) && ...
        y>=handles.display.sliders.s3Pos(2) && ...
        y<handles.display.sliders.s3Pos(2)+handles.display.sliders.s3Pos(4)
    if ~handles.display.sliders.sliderGrab
        if strcmpi(get(handles.display.sliders.s3,'Enable'),'on')
            set(handles.figure.f1,'Pointer','hand')
        else
            set(handles.figure.f1,'Pointer','arrow')
        end
    else
        set(handles.figure.f1,'Pointer','hand')
        turnOffText()
        turnOffCrosshair()
    end
    turnOffText()
    turnOffCrosshair()
    handles.display.cursor.overAxes = false;
elseif (x>=handles.display.sliders.s4Pos(1) && ...
        x<(handles.display.sliders.s4Pos(1)+handles.display.sliders.s4Pos(3)) && ...
        y>=handles.display.sliders.s4Pos(2) && ...
        y<(handles.display.sliders.s4Pos(2)+handles.display.sliders.s4Pos(4))) || ...
        (x>=handles.display.sliders.s5Pos(1) && ...
        x<(handles.display.sliders.s5Pos(1)+handles.display.sliders.s5Pos(3)) && ...
        y>=handles.display.sliders.s5Pos(2) && ...
        y<(handles.display.sliders.s5Pos(2)+handles.display.sliders.s5Pos(4))) || ...
        (x>=handles.display.sliders.s6Pos(1) && ...
        x<(handles.display.sliders.s6Pos(1)+handles.display.sliders.s6Pos(3)) && ...
        y>=handles.display.sliders.s6Pos(2) && ...
        y<(handles.display.sliders.s6Pos(2)+handles.display.sliders.s6Pos(4)))
    set(handles.figure.f1,'Pointer','hand')
    turnOffText()
    turnOffCrosshair()
    handles.display.cursor.overAxes = false;
else
    if ~handles.display.sliders.sliderGrab
        set(handles.figure.f1,'Pointer','arrow')
    else
        set(handles.figure.f1,'Pointer','hand')
    end
    turnOffText()
    turnOffCrosshair()
    handles.display.cursor.overAxes = false;
end

end

function figMotionFcnFused(varargin)

global handles

%ASSUME EVERY UI ELEMENT DEFINED IN NORMALIZED VALUES
units = get(handles.figure.f1,'Units');
set(handles.figure.f1,'Units','normalized');

pt = get(handles.figure.f1,'CurrentPoint');
x = pt(1,1); y = pt(1,2);

set(handles.figure.f1,'Units',units);

if x>=handles.axes.a1Pos(1) && ...
        x<(handles.axes.a1Pos(1)+handles.axes.a1Pos(3)) && ...
        y>=handles.axes.a1Pos(2) && ...
        y<(handles.axes.a1Pos(2)+handles.axes.a1Pos(4))
    if ~handles.display.sliders.sliderGrab
        set(handles.figure.f1,'Pointer','crosshair')
        pt = get(handles.axes.a1,'CurrentPoint');
        pt_x = round(pt(1,1)); pt_y = round(pt(1,2));
        [x,y,z] = getCoordinates(handles.axes.a1,pt_x,pt_y);
        voxelValueFused(x,y,z)
        if handles.display.text.showLocs
            locValue(x,y,z)
        end
        turnOnTextFused()
    else
        set(handles.figure.f1,'Pointer','hand')
        turnOffTextFused()
        turnOffCrosshair()
    end
    moveCrosshair(handles.axes.a1)
    handles.display.cursor.overAxes = true;
elseif x>=handles.axes.a2Pos(1) && ...
        x<(handles.axes.a2Pos(1)+handles.axes.a2Pos(3)) && ...
        y>=handles.axes.a2Pos(2) && ...
        y<(handles.axes.a2Pos(2)+handles.axes.a2Pos(4))
    if ~handles.display.sliders.sliderGrab
        set(handles.figure.f1,'Pointer','crosshair')
        pt = get(handles.axes.a2,'CurrentPoint');
        pt_x = round(pt(1,1)); pt_y = handles.volume.imSz3-round(pt(1,2))+1;
        [x,y,z] = getCoordinates(handles.axes.a2,pt_x,pt_y);
        voxelValueFused(x,y,z)
        if handles.display.text.showLocs
            locValue(x,y,z)
        end
        turnOnTextFused()
    else
        set(handles.figure.f1,'Pointer','hand')
        turnOffTextFused()
        turnOffCrosshair()
    end
    moveCrosshair(handles.axes.a2)
    handles.display.cursor.overAxes = true;
elseif x>=handles.axes.a3Pos(1) && ...
        x<(handles.axes.a3Pos(1)+handles.axes.a3Pos(3)) && ...
        y>=handles.axes.a3Pos(2) && ...
        y<(handles.axes.a3Pos(2)+handles.axes.a3Pos(4))
    if ~handles.display.sliders.sliderGrab
        set(handles.figure.f1,'Pointer','crosshair')
        pt = get(handles.axes.a3,'CurrentPoint');
        pt_x = round(pt(1,1)); pt_y = handles.volume.imSz3-round(pt(1,2))+1;
        [x,y,z] = getCoordinates(handles.axes.a3,pt_x,pt_y);
        voxelValueFused(x,y,z)
        if handles.display.text.showLocs
            locValue(x,y,z)
        end
        turnOnTextFused()
    else
        set(handles.figure.f1,'Pointer','hand')
        turnOffTextFused()
        turnOffCrosshair()
    end
    moveCrosshair(handles.axes.a3)
    handles.display.cursor.overAxes = true;
elseif x>=handles.display.sliders.s1Pos(1) && ...
        x<handles.display.sliders.s1Pos(1)+handles.display.sliders.s1Pos(3) ...
        && y>=handles.display.sliders.s1Pos(2) && ...
        y<handles.display.sliders.s1Pos(2)+handles.display.sliders.s1Pos(4)
    if ~handles.display.sliders.sliderGrab
        if strcmpi(get(handles.display.sliders.s1,'Enable'),'on')
            set(handles.figure.f1,'Pointer','hand')
        else
            set(handles.figure.f1,'Pointer','arrow')
        end
    else
        set(handles.figure.f1,'Pointer','hand')
    end
    turnOffTextFused()
    turnOffCrosshair()
    handles.display.cursor.overAxes = false;
elseif x>=handles.display.sliders.s2Pos(1) && ...
        x<handles.display.sliders.s2Pos(1)+handles.display.sliders.s2Pos(3) && ...
        y>=handles.display.sliders.s2Pos(2) && ...
        y<handles.display.sliders.s2Pos(2)+handles.display.sliders.s2Pos(4)
    if ~handles.display.sliders.sliderGrab
        if strcmpi(get(handles.display.sliders.s2,'Enable'),'on')
            set(handles.figure.f1,'Pointer','hand')
        else
            set(handles.figure.f1,'Pointer','arrow')
        end
    else
        set(handles.figure.f1,'Pointer','hand')
    end
    turnOffTextFused()
    turnOffCrosshair()
    handles.display.cursor.overAxes = false;
elseif x>=handles.display.sliders.s3Pos(1) && ...
        x<handles.display.sliders.s3Pos(1)+handles.display.sliders.s3Pos(3) && ...
        y>=handles.display.sliders.s3Pos(2) && ...
        y<handles.display.sliders.s3Pos(2)+handles.display.sliders.s3Pos(4)
    if ~handles.display.sliders.sliderGrab
        if strcmpi(get(handles.display.sliders.s3,'Enable'),'on')
            set(handles.figure.f1,'Pointer','hand')
        else
            set(handles.figure.f1,'Pointer','arrow')
        end
    else
        set(handles.figure.f1,'Pointer','hand')
        turnOffTextFused()
        turnOffCrosshair()
    end
    turnOffTextFused()
    turnOffCrosshair()
    handles.display.cursor.overAxes = false;
elseif (x>=handles.display.sliders.s4Pos_1(1) && ...
        x<(handles.display.sliders.s4Pos_1(1)+handles.display.sliders.s4Pos_1(3)) && ...
        y>=handles.display.sliders.s4Pos_1(2) && ...
        y<(handles.display.sliders.s4Pos_1(2)+handles.display.sliders.s4Pos_1(4))) || ...
        (x>=handles.display.sliders.s5Pos_1(1) && ...
        x<(handles.display.sliders.s5Pos_1(1)+handles.display.sliders.s5Pos_1(3)) && ...
        y>=handles.display.sliders.s5Pos_1(2) && ...
        y<(handles.display.sliders.s5Pos_1(2)+handles.display.sliders.s5Pos_1(4))) || ...
        (x>=handles.display.sliders.s4Pos_2(1) && ...
        x<(handles.display.sliders.s4Pos_2(1)+handles.display.sliders.s4Pos_2(3)) && ...
        y>=handles.display.sliders.s4Pos_2(2) && ...
        y<(handles.display.sliders.s4Pos_2(2)+handles.display.sliders.s4Pos_2(4))) || ...
        (x>=handles.display.sliders.s5Pos_2(1) && ...
        x<(handles.display.sliders.s5Pos_2(1)+handles.display.sliders.s5Pos_2(3)) && ...
        y>=handles.display.sliders.s5Pos_2(2) && ...
        y<(handles.display.sliders.s5Pos_2(2)+handles.display.sliders.s5Pos_2(4))) || ...
        (x>=handles.display.sliders.s6Pos_1(1) && ...
        x<(handles.display.sliders.s6Pos_1(1)+handles.display.sliders.s6Pos_1(3)) && ...
        y>=handles.display.sliders.s6Pos_1(2) && ...
        y<(handles.display.sliders.s6Pos_1(2)+handles.display.sliders.s6Pos_1(4)))
    set(handles.figure.f1,'Pointer','hand')
    turnOffTextFused()
    turnOffCrosshair()
    handles.display.cursor.overAxes = false;
else
    if ~handles.display.sliders.sliderGrab
        set(handles.figure.f1,'Pointer','arrow')
    else
        set(handles.figure.f1,'Pointer','hand')
    end
    turnOffTextFused()
    turnOffCrosshair()
    handles.display.cursor.overAxes = false;
end

end

function [x,y,z] = getCoordinates(a,x,y)

global handles

imSzX = handles.volume.imSz1;
imSzY = handles.volume.imSz2;
imSzZ = handles.volume.imSz3;

switch a
    case handles.axes.a1
        a1xlim = get(handles.axes.a1,'XLim');
        a1ylim = get(handles.axes.a1,'YLim');
        if x<1 || x<a1xlim(1)
            x = max(1,a1xlim(1));
        elseif x>imSzX || x>a1xlim(2)
            x = min(imSzX,a1xlim(2));
        end
        if y<1 || y<a1ylim(1)
            y = max(1,a1ylim(1));
        elseif y>imSzY || y>a1ylim(2)
            y = min(imSzY,a1ylim(2));
        end
        %         if x<1 || x>imSzX
        %             x = handles.display.images.x;
        %         end
        %         if y<1 || y>imSzY
        %             y = imSzY-handles.display.images.y+1;
        %         end
        x = round(x);
        y = round(imSzY-y+1);
        z = round(handles.display.images.z);
    case handles.axes.a2
        y = handles.volume.imSz3-y+1;
        a2xlim = get(handles.axes.a2,'XLim');
        a2ylim = get(handles.axes.a2,'YLim');
        if x<1 || x<a2xlim(1)
            x = max(1,a2xlim(1));
        elseif x>imSzX || x>a2xlim(2)
            x = min(imSzX,a2xlim(2));
        end
        if y<1 || y<a2ylim(1)
            y = max(1,a2ylim(1));
        elseif y>imSzZ || y>a2ylim(2)
            y = min(imSzZ,a2ylim(2));
        end
        %         if x<1 || x>imSzX
        %             x = handles.display.images.x;
        %         end
        %         if y<1 || y>imSzZ
        %             y = imSzZ-handles.display.images.z+1;
        %         end
        x = round(x);
        z = round(y);
        y = handles.display.images.y;
    case handles.axes.a3
        y = handles.volume.imSz3-y+1;
        a3xlim = get(handles.axes.a3,'XLim');
        a3ylim = get(handles.axes.a3,'YLim');
        if x<1 || x<a3xlim(1)
            x = max(1,a3xlim(1));
        elseif x>imSzY || x>a3xlim(2)
            x = min(imSzY,a3xlim(2));
        end
        if y<1 || y<a3ylim(1)
            y = max(1,a3ylim(1));
        elseif y>imSzZ || y>a3ylim(2)
            y = min(imSzZ,a3ylim(2));
        end
        %         if x<1 || x>imSzY
        %             x = handles.display.images.y;
        %         end
        %         if y<1 || y>imSzZ
        %             y = imSzZ-handles.display.images.z+1;
        %         end
        z = round(y);
        y = round(x);
        x = handles.display.images.x;
end

if x<handles.display.images.xLimits(1), x = handles.display.images.xLimits(1); end
if x>handles.display.images.xLimits(2), x = handles.display.images.xLimits(2); end
if y<handles.display.images.yLimits(1), y = handles.display.images.yLimits(1); end
if y>handles.display.images.yLimits(2), y = handles.display.images.yLimits(2); end
if z<handles.display.images.zLimits(1), z = handles.display.images.zLimits(1); end
if z>handles.display.images.zLimits(2), z = handles.display.images.zLimits(2); end

end

function voxelValue(x,y,z)

global handles

str = formatValueString(handles.volume.vol(x,y,z)/handles.display.images.sclFctr);

set(handles.display.text.vTxt,'String',str)

end

function voxelValueFused(x,y,z)

global handles

str1 = formatValueString(handles.volume.vol_1(x,y,z)/handles.display.images.sclFctr_1);
str2 = formatValueString(handles.volume.vol_2(x,y,z)/handles.display.images.sclFctr_2);

set(handles.display.text.vTxt_1,'String',str1)
set(handles.display.text.vTxt_2,'String',str2)

end

function coordinatesValue(x,y,z)

global handles

set(handles.display.text.cTxt,'String',sprintf('%0.0f, %0.0f, %0.0f',x,y,z))

end

function locValue(x,y,z)

global handles

set(handles.display.text.locTxtX,'String',sprintf('%0.2f mm',handles.volume.locX(x,y,z)))
set(handles.display.text.locTxtY,'String',sprintf('%0.2f mm',handles.volume.locY(x,y,z)))
set(handles.display.text.locTxtZ,'String',sprintf('%0.2f mm',handles.volume.locZ(x,y,z)))

end

function toggleFcn(varargin)

global handles

if strcmpi(get(gcbf,'SelectionType'),'open')
    if ~handles.display.images.fusedImages && ~strcmpi(handles.volume.modality,'CT')
        invertCM()
    end
    % elseif strcmpi(get(gcbf,'SelectionType'),'alt')
    %     toggleCrosshair()
elseif strcmpi(get(gcbf,'SelectionType'),'normal')
    
    %THIS ONLY NEEDED IF CURSOR OVER (NON-HITTABLE) ARROWHEADS AT IMAGE EDGES
    units = get(handles.figure.f1,'Units');
    
    %ASSUME EVERY UI ELEMENT DEFINED IN NORMALIZED VALUES
    set(handles.figure.f1,'Units','normalized');
    pt = get(handles.figure.f1,'CurrentPoint');
    x = pt(1,1); y = pt(1,2);
    set(handles.figure.f1,'Units',units);
    
    if x>=handles.axes.a1Pos(1) && ...
            x<(handles.axes.a1Pos(1)+handles.axes.a1Pos(3)) && ...
            y>=handles.axes.a1Pos(2) && ...
            y<(handles.axes.a1Pos(2)+handles.axes.a1Pos(4))
        a1ButtonDownFcn(handles.axes.a1)
    elseif x>=handles.axes.a2Pos(1) && ...
            x<(handles.axes.a2Pos(1)+handles.axes.a2Pos(3)) && ...
            y>=handles.axes.a2Pos(2) && ...
            y<(handles.axes.a2Pos(2)+handles.axes.a2Pos(4))
        a2ButtonDownFcn(handles.axes.a2)
    elseif x>=handles.axes.a3Pos(1) && ...
            x<(handles.axes.a3Pos(1)+handles.axes.a3Pos(3)) && ...
            y>=handles.axes.a3Pos(2) && ...
            y<(handles.axes.a3Pos(2)+handles.axes.a3Pos(4))
        a3ButtonDownFcn(handles.axes.a3)
    end
    
end

end

function moveCrosshair(a)

global handles

XLim = xlim(handles.axes.a1);
YLim = ylim(handles.axes.a1);
ZLim = ylim(handles.axes.a2);
%ACCOUNT FOR Y FLIP
try
    a3YLim = flip(handles.volume.imSz2+1-YLim);
catch
    a3YLim = fliplr(handles.volume.imSz2+1-YLim);
end

switch a
    case handles.axes.a1
        pt = get(handles.axes.a1,'CurrentPoint');
        pt_x = round(pt(1,1)); pt_y = round(pt(1,2));
        %         if (pt_x>=1 && pt_x<=imSzX) || (pt_y>=1 && pt_y<=imSzY)
        %a1
        [x,y,~] = getCoordinates(a,pt_x,pt_y);
        crossX = x;
        crossY = handles.volume.imSz2-y+1;
        xGap = round(handles.display.crosshairs.crossGap/handles.volume.voxDimX/2);
        yGap = round(handles.display.crosshairs.crossGap/handles.volume.voxDimY/2);
        set(handles.display.crosshairs.a1crossH1,'XData',(XLim(1)-1:x-xGap),'YData',crossY*(ones(size((XLim(1)-1:x-xGap)))))
        set(handles.display.crosshairs.a1crossH2,'XData',(crossX+xGap:XLim(2)+1),'YData',crossY*(ones(size((crossX+xGap:XLim(2)+1)))))
        set(handles.display.crosshairs.a1crossV1,'XData',crossX*(ones(size((YLim(1)-1:crossY-yGap)))),'YData',(YLim(1)-1:crossY-yGap))
        set(handles.display.crosshairs.a1crossV2,'XData',crossX*(ones(size((crossY+yGap:YLim(2)+1)))),'YData',(crossY+yGap:YLim(2)+1))
        %a2
        crossY = handles.display.images.z;
        yGap = round(handles.display.crosshairs.crossGap/handles.volume.voxDimZ/2);
        set(handles.display.crosshairs.a2crossH1,'XData',(XLim(1)-1:x-xGap),'YData',crossY*(ones(size((XLim(1)-1:x-xGap)))))
        set(handles.display.crosshairs.a2crossH2,'XData',(crossX+xGap:XLim(2)+1),'YData',crossY*(ones(size((crossX+xGap:XLim(2)+1)))))
        set(handles.display.crosshairs.a2crossV1,'XData',crossX*ones(size(ZLim(1)-1:crossY-yGap)),'YData',(ZLim(1)-1:crossY-yGap))
        set(handles.display.crosshairs.a2crossV2,'XData',crossX*ones(size(crossY+yGap:ZLim(2)+1)),'YData',(crossY+yGap:ZLim(2)+1))
        %a3
        crossX = y;
        xGap = round(handles.display.crosshairs.crossGap/handles.volume.voxDimY/2);
        set(handles.display.crosshairs.a3crossH1,'XData',(a3YLim(1)-1:crossX-xGap),'YData',crossY*(ones(size((a3YLim(1)-1:crossX-xGap)))))
        set(handles.display.crosshairs.a3crossH2,'XData',(crossX+xGap:a3YLim(2)+1),'YData',crossY*(ones(size((crossX+xGap:a3YLim(2)+1)))))
        set(handles.display.crosshairs.a3crossV1,'XData',crossX*ones(size(ZLim(1)-1:crossY-yGap)),'YData',(ZLim(1)-1:crossY-yGap))
        set(handles.display.crosshairs.a3crossV2,'XData',crossX*ones(size(crossY+yGap:ZLim(2)+1)),'YData',(crossY+yGap:ZLim(2)+1))
        %         end
    case handles.axes.a2
        pt = get(handles.axes.a2,'CurrentPoint');
        pt_x = round(pt(1,1)); pt_y = round(pt(1,2));
        %         if (pt_x>=1 && pt_x<=imSzX) || (pt_y>=1 && pt_y<=imSzZ)
        %a2
        [x,~,y] = getCoordinates(a,pt_x,handles.volume.imSz3-pt_y+1);
        crossX = x;
        crossY = y;
        xGap = round(handles.display.crosshairs.crossGap/handles.volume.voxDimX/2);
        yGap = round(handles.display.crosshairs.crossGap/handles.volume.voxDimZ/2);
        set(handles.display.crosshairs.a2crossH1,'XData',(XLim(1)-1:x-xGap),'YData',crossY*(ones(size((XLim(1)-1:x-xGap)))))
        set(handles.display.crosshairs.a2crossH2,'XData',(crossX+xGap:XLim(2)+1),'YData',crossY*(ones(size((crossX+xGap:XLim(2)+1)))))
        set(handles.display.crosshairs.a2crossV1,'XData',crossX*ones(size(ZLim(1)-1:crossY-yGap)),'YData',(ZLim(1)-1:crossY-yGap))
        set(handles.display.crosshairs.a2crossV2,'XData',crossX*ones(size(crossY+yGap:ZLim(2)+1)),'YData',(crossY+yGap:ZLim(2)+1))
        %a3
        crossX = handles.display.images.y;
        xGap = round(handles.display.crosshairs.crossGap/handles.volume.voxDimY/2);
        set(handles.display.crosshairs.a3crossH1,'XData',(a3YLim(1)-1:crossX-xGap),'YData',crossY*(ones(size((a3YLim(1)-1:crossX-xGap)))))
        set(handles.display.crosshairs.a3crossH2,'XData',(crossX+xGap:a3YLim(2)+1),'YData',crossY*(ones(size((crossX+xGap:a3YLim(2)+1)))))
        set(handles.display.crosshairs.a3crossV1,'XData',crossX*ones(size(ZLim(1)-1:crossY-yGap)),'YData',(ZLim(1)-1:crossY-yGap))
        set(handles.display.crosshairs.a3crossV2,'XData',crossX*ones(size(crossY+yGap:ZLim(2)+1)),'YData',(crossY+yGap:ZLim(2)+1))
        %a1
        crossX = x;
        crossY = handles.volume.imSz2-handles.display.images.y+1;
        xGap = round(handles.display.crosshairs.crossGap/handles.volume.voxDimX/2);
        yGap = round(handles.display.crosshairs.crossGap/handles.volume.voxDimY/2);
        set(handles.display.crosshairs.a1crossH1,'XData',(XLim(1)-1:x-xGap),'YData',crossY*(ones(size((XLim(1)-1:x-xGap)))))
        set(handles.display.crosshairs.a1crossH2,'XData',(crossX+xGap:XLim(2)+1),'YData',crossY*(ones(size((crossX+xGap:XLim(2)+1)))))
        set(handles.display.crosshairs.a1crossV1,'XData',crossX*(ones(size((YLim(1)-1:crossY-yGap)))),'YData',(YLim(1)-1:crossY-yGap))
        set(handles.display.crosshairs.a1crossV2,'XData',crossX*(ones(size((crossY+yGap:YLim(2)+1)))),'YData',(crossY+yGap:YLim(2)+1))
        %         end
    case handles.axes.a3
        pt = get(handles.axes.a3,'CurrentPoint');
        pt_x = round(pt(1,1)); pt_y = round(pt(1,2));
        %         if (pt_x>=1 && pt_x<=imSzY) || (pt_y>=1 && pt_y<=imSzZ)
        %a3
        [~,y,z] = getCoordinates(a,pt_x,handles.volume.imSz3-pt_y+1);
        crossX = y;
        crossY = z;
        xGap = round(handles.display.crosshairs.crossGap/handles.volume.voxDimY/2);
        yGap = round(handles.display.crosshairs.crossGap/handles.volume.voxDimZ/2);
        set(handles.display.crosshairs.a3crossH1,'XData',(a3YLim(1)-1:crossX-xGap),'YData',crossY*(ones(size((a3YLim(1)-1:crossX-xGap)))))
        set(handles.display.crosshairs.a3crossH2,'XData',(crossX+xGap:a3YLim(2)+1),'YData',crossY*(ones(size((crossX+xGap:a3YLim(2)+1)))))
        set(handles.display.crosshairs.a3crossV1,'XData',crossX*ones(size(ZLim(1)-1:crossY-yGap)),'YData',(ZLim(1)-1:crossY-yGap))
        set(handles.display.crosshairs.a3crossV2,'XData',crossX*ones(size(crossY+yGap:ZLim(2)+1)),'YData',(crossY+yGap:ZLim(2)+1))
        %a2
        crossX = handles.display.images.x;
        xGap = round(handles.display.crosshairs.crossGap/handles.volume.voxDimX/2);
        set(handles.display.crosshairs.a2crossH1,'XData',(XLim(1)-1:crossX-xGap),'YData',crossY*(ones(size((XLim(1)-1:crossX-xGap)))))
        set(handles.display.crosshairs.a2crossH2,'XData',(crossX+xGap:XLim(2)+1),'YData',crossY*(ones(size((crossX+xGap:XLim(2)+1)))))
        set(handles.display.crosshairs.a2crossV1,'XData',crossX*ones(size(ZLim(1)-1:crossY-yGap)),'YData',(ZLim(1)-1:crossY-yGap))
        set(handles.display.crosshairs.a2crossV2,'XData',crossX*ones(size(crossY+yGap:ZLim(2)+1)),'YData',(crossY+yGap:ZLim(2)+1))
        %a1
        crossY = handles.volume.imSz2-y+1;
        yGap = round(handles.display.crosshairs.crossGap/handles.volume.voxDimY/2);
        set(handles.display.crosshairs.a1crossH1,'XData',(XLim(1)-1:crossX-xGap),'YData',crossY*(ones(size((XLim(1)-1:crossX-xGap)))))
        set(handles.display.crosshairs.a1crossH2,'XData',(crossX+xGap:XLim(2)+1),'YData',crossY*(ones(size((crossX+xGap:XLim(2)+1)))))
        set(handles.display.crosshairs.a1crossV1,'XData',crossX*(ones(size((YLim(1)-1:crossY-yGap)))),'YData',(YLim(1)-1:crossY-yGap))
        set(handles.display.crosshairs.a1crossV2,'XData',crossX*(ones(size((crossY+yGap:YLim(2)+1)))),'YData',(crossY+yGap:YLim(2)+1))
        %         end
end

if handles.display.crosshairs.crossBoo && ~handles.display.sliders.sliderGrab
    set(handles.display.crosshairs.a1crossH1,'Visible','on')
    set(handles.display.crosshairs.a1crossH2,'Visible','on')
    set(handles.display.crosshairs.a1crossV1,'Visible','on')
    set(handles.display.crosshairs.a1crossV2,'Visible','on')
    set(handles.display.crosshairs.a2crossH1,'Visible','on')
    set(handles.display.crosshairs.a2crossH2,'Visible','on')
    set(handles.display.crosshairs.a2crossV1,'Visible','on')
    set(handles.display.crosshairs.a2crossV2,'Visible','on')
    set(handles.display.crosshairs.a3crossH1,'Visible','on')
    set(handles.display.crosshairs.a3crossH2,'Visible','on')
    set(handles.display.crosshairs.a3crossV1,'Visible','on')
    set(handles.display.crosshairs.a3crossV2,'Visible','on')
end

end

function toggleCrosshair

global handles

if handles.display.crosshairs.crossBoo==0
    handles.display.crosshairs.crossBoo = true;
    %do not show crosshairs if mouse cursor is outside of axes
    pt = get(handles.figure.f1,'CurrentPoint');
    x = pt(1,1); y = pt(1,2);
    
    a1Pos = get(handles.axes.a1,'Position');
    a2Pos = get(handles.axes.a2,'Position');
    a3Pos = get(handles.axes.a3,'Position');
    
    if (x>=a1Pos(1) && x<(a1Pos(1)+a1Pos(3)) && y>=a1Pos(2)+1 && y<(a1Pos(2)+a1Pos(4))) || ...
            (x>=a2Pos(1) && x<(a2Pos(1)+a2Pos(3)) && y>=a2Pos(2) && y<(a2Pos(2)+a2Pos(4))) || ...
            (x>=a3Pos(1) && x<(a3Pos(1)+a3Pos(3)) && y>=a3Pos(2) && y<(a3Pos(2)+a3Pos(4)))
        turnOnCrosshair()
    end
else
    handles.display.crosshairs.crossBoo = false;
    turnOffCrosshair()
end

end

function turnOnCrosshair

global handles

set(handles.display.crosshairs.a1crossH1,'Visible','on')
set(handles.display.crosshairs.a1crossH2,'Visible','on')
set(handles.display.crosshairs.a1crossV1,'Visible','on')
set(handles.display.crosshairs.a1crossV2,'Visible','on')
set(handles.display.crosshairs.a2crossH1,'Visible','on')
set(handles.display.crosshairs.a2crossH2,'Visible','on')
set(handles.display.crosshairs.a2crossV1,'Visible','on')
set(handles.display.crosshairs.a2crossV2,'Visible','on')
set(handles.display.crosshairs.a3crossH1,'Visible','on')
set(handles.display.crosshairs.a3crossH2,'Visible','on')
set(handles.display.crosshairs.a3crossV1,'Visible','on')
set(handles.display.crosshairs.a3crossV2,'Visible','on')

end

function turnOffCrosshair

global handles

set(handles.display.crosshairs.a1crossH1,'Visible','off')
set(handles.display.crosshairs.a1crossH2,'Visible','off')
set(handles.display.crosshairs.a1crossV1,'Visible','off')
set(handles.display.crosshairs.a1crossV2,'Visible','off')
set(handles.display.crosshairs.a2crossH1,'Visible','off')
set(handles.display.crosshairs.a2crossH2,'Visible','off')
set(handles.display.crosshairs.a2crossV1,'Visible','off')
set(handles.display.crosshairs.a2crossV2,'Visible','off')
set(handles.display.crosshairs.a3crossH1,'Visible','off')
set(handles.display.crosshairs.a3crossH2,'Visible','off')
set(handles.display.crosshairs.a3crossV1,'Visible','off')
set(handles.display.crosshairs.a3crossV2,'Visible','off')

end

function moveArrows

global handles

a1hidx = handles.display.arrows.a1HarrowLocs(handles.volume.imSz2-handles.display.images.y+1);
a1vidx = handles.display.arrows.a1VarrowLocs(handles.display.images.x);
a2hidx = handles.display.arrows.a2HarrowLocs(handles.display.images.z);
a2vidx = handles.display.arrows.a2VarrowLocs(handles.display.images.x);
a3hidx = handles.display.arrows.a3HarrowLocs(handles.display.images.z);
a3vidx = handles.display.arrows.a3VarrowLocs(handles.display.images.y);

set(handles.display.arrows.a1arrowH1,'Y',[a1hidx a1hidx])
set(handles.display.arrows.a1arrowH2,'Y',[a1hidx a1hidx])
set(handles.display.arrows.a1arrowV1,'X',[a1vidx a1vidx])
set(handles.display.arrows.a1arrowV2,'X',[a1vidx a1vidx])
set(handles.display.arrows.a2arrowH1,'Y',[a2hidx a2hidx])
set(handles.display.arrows.a2arrowH2,'Y',[a2hidx a2hidx])
set(handles.display.arrows.a2arrowV1,'X',[a2vidx a2vidx])
set(handles.display.arrows.a2arrowV2,'X',[a2vidx a2vidx])
set(handles.display.arrows.a3arrowH1,'Y',[a3hidx a3hidx])
set(handles.display.arrows.a3arrowH2,'Y',[a3hidx a3hidx])
set(handles.display.arrows.a3arrowV1,'X',[a3vidx a3vidx])
set(handles.display.arrows.a3arrowV2,'X',[a3vidx a3vidx])

end

function turnOnArrows

global handles

set(handles.display.arrows.a1arrowH1,'Visible','on')
set(handles.display.arrows.a1arrowH2,'Visible','on')
set(handles.display.arrows.a1arrowV1,'Visible','on')
set(handles.display.arrows.a1arrowV2,'Visible','on')
set(handles.display.arrows.a2arrowH1,'Visible','on')
set(handles.display.arrows.a2arrowH2,'Visible','on')
set(handles.display.arrows.a2arrowV1,'Visible','on')
set(handles.display.arrows.a2arrowV2,'Visible','on')
set(handles.display.arrows.a3arrowH1,'Visible','on')
set(handles.display.arrows.a3arrowH2,'Visible','on')
set(handles.display.arrows.a3arrowV1,'Visible','on')
set(handles.display.arrows.a3arrowV2,'Visible','on')

end

function turnOffArrows

global handles

set(handles.display.arrows.a1arrowH1,'Visible','off')
set(handles.display.arrows.a1arrowH2,'Visible','off')
set(handles.display.arrows.a1arrowV1,'Visible','off')
set(handles.display.arrows.a1arrowV2,'Visible','off')
set(handles.display.arrows.a2arrowH1,'Visible','off')
set(handles.display.arrows.a2arrowH2,'Visible','off')
set(handles.display.arrows.a2arrowV1,'Visible','off')
set(handles.display.arrows.a2arrowV2,'Visible','off')
set(handles.display.arrows.a3arrowH1,'Visible','off')
set(handles.display.arrows.a3arrowH2,'Visible','off')
set(handles.display.arrows.a3arrowV1,'Visible','off')
set(handles.display.arrows.a3arrowV2,'Visible','off')

end

function turnOnText

global handles

set(handles.display.text.vTxt,'Visible','on')
set(handles.display.text.cTxt,'Visible','on')
if handles.display.text.showLocs
    set(handles.display.text.locTxtX,'Visible','on')
    set(handles.display.text.locTxtY,'Visible','on')
    set(handles.display.text.locTxtZ,'Visible','on')
end

end

function turnOffText

global handles

set(handles.display.text.vTxt,'Visible','off')
set(handles.display.text.cTxt,'Visible','off')
if handles.display.text.showLocs
    set(handles.display.text.locTxtX,'Visible','off')
    set(handles.display.text.locTxtY,'Visible','off')
    set(handles.display.text.locTxtZ,'Visible','off')
end

end

function turnOnTextFused

global handles

set(handles.display.text.vTxt_1,'Visible','on')
set(handles.display.text.vTxt_2,'Visible','on')
if handles.display.text.showLocs
    set(handles.display.text.locTxtX,'Visible','on')
    set(handles.display.text.locTxtY,'Visible','on')
    set(handles.display.text.locTxtZ,'Visible','on')
end

end

function turnOffTextFused

global handles

set(handles.display.text.vTxt_1,'Visible','off')
set(handles.display.text.vTxt_2,'Visible','off')
if handles.display.text.showLocs
    set(handles.display.text.locTxtX,'Visible','off')
    set(handles.display.text.locTxtY,'Visible','off')
    set(handles.display.text.locTxtZ,'Visible','off')
end

end

function labelPatientOrientation

global handles

figSz = get(handles.figure.f1,'Position');
figSz = [figSz(3) figSz(4)];

txtCol = [0.9290 0.6940 0.1250];
txtSz = 0.006*figSz(1)*handles.figure.DPIscl;
if handles.figure.DPIscl~=1 %high DPI scaling
    txtSz = txtSz * 0.7;
end

rot = [handles.volume.IOP(1:3) handles.volume.IOP(4:6) ...
    cross(handles.volume.IOP(1:3),handles.volume.IOP(4:6))]';

%make a matrix of orientation labels to match patient rotation
%the higher the matrix dimensionality, the stricter the label criteria
labMatDim = 51; %odd number
labMatDim = labMatDim-mod(labMatDim,2)+1; %assure odd number
labMat = cell(labMatDim,labMatDim,labMatDim);
idx1 = 1; idx2 = round(size(labMat,1)/2); idx3 = size(labMat,1);
%a numerical matrix to track positions for label positioning
labPtrMat = zeros(labMatDim,labMatDim,labMatDim);
%and pointer matrix, convention: a1t1,a1t2,a2t1,a2t2,a3t1,a3t2 -> 1,2,1,3,4,3
ptr = [...
    idx2 idx2 idx1;...
    idx2 idx2 idx3;...
    idx2 idx1 idx2;...
    idx2 idx3 idx2;...
    idx1 idx2 idx2;...
    idx3 idx2 idx2];
%and rotate
ptr = round(rot * (ptr' - repmat(idx2,3,6)) + repmat(idx2,3,6));
%new orientation
labMat(ptr(1,1),ptr(2,1),ptr(3,1)) = {' F '};
labMat(ptr(1,2),ptr(2,2),ptr(3,2)) = {' H '};
labMat(ptr(1,3),ptr(2,3),ptr(3,3)) = {' A '};
labMat(ptr(1,4),ptr(2,4),ptr(3,4)) = {' P '};
labMat(ptr(1,5),ptr(2,5),ptr(3,5)) = {' R '};
labMat(ptr(1,6),ptr(2,6),ptr(3,6)) = {' L '};
a1t1 = labMat{idx1,idx2,idx2};
a2t1 = labMat{idx1,idx2,idx2};
a1t2 = labMat{idx2,idx3,idx2};
a3t1 = labMat{idx2,idx1,idx2};
a2t2 = labMat{idx2,idx2,idx1};
a3t2 = labMat{idx2,idx2,idx1};
%make pointer indeces odd so we know if sum is even, it was overlapped
labPtrMat(ptr(1,1),ptr(2,1),ptr(3,1)) = 1;
labPtrMat(ptr(1,2),ptr(2,2),ptr(3,2)) = 3;
labPtrMat(ptr(1,3),ptr(2,3),ptr(3,3)) = 5;
labPtrMat(ptr(1,4),ptr(2,4),ptr(3,4)) = 7;
labPtrMat(ptr(1,5),ptr(2,5),ptr(3,5)) = 9;
labPtrMat(ptr(1,6),ptr(2,6),ptr(3,6)) = 11;

%check relative direction of bed ranges
%X
if handles.volume.bedRange1(1)>handles.volume.bedRange1(2)
    try
        labPtrMat = flip(labPtrMat,1);
    catch
        labPtrMat = flipdim(labPtrMat,1);
    end
end
%Y
if handles.volume.bedRange2(1)>handles.volume.bedRange2(2)
    try
        labPtrMat = flip(labPtrMat,2);
    catch
        labPtrMat = flipdim(labPtrMat,2);
    end
end
%Z
if handles.volume.bedRange3(1)>handles.volume.bedRange3(2)
    try
        labPtrMat = flip(labPtrMat,3);
    catch
        labPtrMat = flipdim(labPtrMat,3);
    end
end

set(handles.axes.a1,'Units','points')
set(handles.axes.a2,'Units','points')
set(handles.axes.a3,'Units','points')
ax1Pos = get(handles.axes.a1,'Position');
ax2Pos = get(handles.axes.a2,'Position');
ax3Pos = get(handles.axes.a3,'Position');
set(handles.axes.a1,'Units','pixels')
set(handles.axes.a2,'Units','pixels')
set(handles.axes.a3,'Units','pixels')
%offsets for labels, these may be changed if using rotated frame of reference
a1t1Offset = 0;
a1t2Offset = 0;
a2t1Offset = 0;
a2t2Offset = 0;
a3t1Offset = 0;
a3t2Offset = 0;
axHorizLabPad = 2*txtSz+9; %axes edge padding
axVertLabPad = 2*txtSz+11; %axes edge padding

%sum in each direction
sum1 = squeeze(sum(labPtrMat,1));
sum2 = squeeze(sum(labPtrMat,2));
sum3 = squeeze(sum(labPtrMat,3));
%the first 2 found indeces in each direction can be used for orientation
%axes 1
sumsum = sum(sum3,2); %sum rows
temp = find(sumsum~=0);
if numel(find(sum3~=0))==5 && numel(temp)<=4 && all(mod(sumsum,2)==0) %45 degree rotation
    temp = find(sum3~=0);
    a1t1 = getLabelFromIndex(sum3(temp(end-1)));
    a1t2 = getLabelFromIndex(sum3(temp(end)));
    %offsets
    idx = find(sum3==sum3(temp(end-1)));
    x = abs(mod(idx,labMatDim)-round(labMatDim/2));
    y = -(ceil(idx/labMatDim)-round(labMatDim/2));
    a1t1Offset= y/x*(ax1Pos(4)-axVertLabPad)/2;
    idx = find(sum3==sum3(temp(end)));
    x = mod(idx,labMatDim)-round(labMatDim/2);
    y = abs(ceil(idx/labMatDim)-round(labMatDim/2));
    a1t2Offset = x/y*(ax1Pos(3)-axHorizLabPad)/2;
else
    if mod(sumsum(temp(1)),2)==1 && numel(temp)>2
        a1t1 = getLabelFromIndex(sumsum(temp(1)));
        %offset
        idx = find(sum3==sumsum(temp(1)));
        x = abs(mod(idx,labMatDim)-round(labMatDim/2));
        y = -(ceil(idx/labMatDim)-round(labMatDim/2));
        a1t1Offset= y/x*(ax1Pos(4)-axVertLabPad)/2;
    end
    sumsum = sum(sum3,1); %sum columns
    temp = find(sumsum~=0);
    if mod(sumsum(temp(end)),2)==1 && numel(temp)>2
        a1t2 = getLabelFromIndex(sumsum(temp(end)));
        %offset
        idx = find(sum3==sumsum(temp(end)));
        x = mod(idx,labMatDim)-round(labMatDim/2);
        y = abs(ceil(idx/labMatDim)-round(labMatDim/2));
        a1t2Offset = x/y*(ax1Pos(3)-axHorizLabPad)/2;
    end
end
%axes 2
sumsum = sum(sum2,2); %sum rows
temp = find(sumsum~=0);
if numel(find(sum2~=0))==5 && numel(temp)<=4 && all(mod(sumsum,2)==0) %45 degree rotation
    temp = find(sum2~=0);
    a2t1 = getLabelFromIndex(sum2(temp(1)));
    a2t2 = getLabelFromIndex(sum2(temp(2)));
    %offsets
    idx = find(sum2==sum2(temp(1)));
    x = abs(mod(idx,labMatDim)-round(labMatDim/2));
    y = ceil(idx/labMatDim)-round(labMatDim/2);
    a2t1Offset= y/x*(ax2Pos(4)-axVertLabPad)/2;
    idx = find(sum2==sum2(temp(2)));
    x = mod(idx,labMatDim)-round(labMatDim/2);
    y = abs(ceil(idx/labMatDim)-round(labMatDim/2));
    a2t2Offset = x/y*(ax2Pos(3)-axHorizLabPad)/2;
else
    if mod(sumsum(temp(1)),2)==1 && numel(temp)>2
        a2t1 = getLabelFromIndex(sumsum(temp(1)));
        %offset
        idx = find(sum2==sumsum(temp(1)));
        x = abs(mod(idx,labMatDim)-round(labMatDim/2));
        y = ceil(idx/labMatDim)-round(labMatDim/2);
        a2t1Offset= y/x*(ax2Pos(4)-axVertLabPad)/2;
    end
    sumsum = sum(sum2,1); %sum columns
    temp = find(sumsum~=0);
    if mod(sumsum(temp(1)),2)==1 && numel(temp)>2
        a2t2 = getLabelFromIndex(sumsum(temp(1)));
        %offset
        idx = find(sum2==sumsum(temp(1)));
        x = mod(idx,labMatDim)-round(labMatDim/2);
        y = abs(ceil(idx/labMatDim)-round(labMatDim/2));
        a2t2Offset = x/y*(ax2Pos(3)-axHorizLabPad)/2;
    end
end
%axes 3
sumsum = sum(sum1,2); %sum rows
temp = find(sumsum~=0);
if numel(find(sum1~=0))==5 && numel(temp)<=4 && all(mod(sumsum,2)==0) %45 degree rotation
    temp = find(sum1~=0);
    a3t1 = getLabelFromIndex(sum1(temp(1)));
    a3t2 = getLabelFromIndex(sum1(temp(2)));
    %offsets
    idx = find(sum1==sum1(temp(1)));
    x = abs(mod(idx,labMatDim)-round(labMatDim/2));
    y = ceil(idx/labMatDim)-round(labMatDim/2);
    a3t1Offset= y/x*(ax3Pos(4)-axVertLabPad)/2;
    idx = find(sum1==sum1(temp(2)));
    x = mod(idx,labMatDim)-round(labMatDim/2);
    y = abs(ceil(idx/labMatDim)-round(labMatDim/2));
    a3t2Offset = x/y*(ax3Pos(3)-axHorizLabPad)/2;
else
    if mod(sumsum(temp(1)),2)==1 && numel(temp)>2
        a3t1 = getLabelFromIndex(sumsum(temp(1)));
        %offset
        idx = find(sum1==sumsum(temp(1)));
        x = abs(mod(idx,labMatDim)-round(labMatDim/2));
        y = ceil(idx/labMatDim)-round(labMatDim/2);
        a3t1Offset= y/x*(ax3Pos(4)-axVertLabPad)/2;
    end
    sumsum = sum(sum1,1); %sum columns
    temp = find(sumsum~=0);
    if mod(sumsum(temp(1)),2)==1 && numel(temp)>2
        a3t2 = getLabelFromIndex(sumsum(temp(1)));
        %offset
        idx = find(sum1==sumsum(temp(1)));
        x = mod(idx,labMatDim)-round(labMatDim/2);
        y = abs(ceil(idx/labMatDim)-round(labMatDim/2));
        a3t2Offset = x/y*(ax3Pos(3)-axHorizLabPad)/2;
    end
end

%axes 1 title
if a2t2Offset==0 && a3t2Offset==0
    if ismember('H',a2t2) || ismember('F',a2t2)
        handles.axes.a1title = {'Transaxial';''};
    elseif ismember('A',a2t2) || ismember('P',a2t2)
        handles.axes.a1title = {'Coronal';''};
    elseif ismember('L',a2t2) || ismember('R',a2t2)
        handles.axes.a1title = {'Sagittal';''};
    end
end
%axes 2 title
if a1t2Offset==0 && a3t1Offset==0
    if ismember('H',a1t2) || ismember('F',a1t2)
        handles.axes.a2title = {'Transaxial';''};
    elseif ismember('A',a1t2) || ismember('P',a1t2)
        handles.axes.a2title = {'Coronal';''};
    elseif ismember('L',a1t2) || ismember('R',a1t2)
        handles.axes.a2title = {'Sagittal';''};
    end
end
%axes 3 title
if a1t1Offset==0 && a2t1Offset==0
    if ismember('H',a1t1) || ismember('F',a1t1)
        handles.axes.a3title = {'Transaxial';''};
    elseif ismember('A',a1t1) || ismember('P',a1t1)
        handles.axes.a3title = {'Coronal';''};
    elseif ismember('L',a1t1) || ismember('R',a1t1)
        handles.axes.a3title = {'Sagittal';''};
    end
end

%label positions
a1t1Pos = [8 ax1Pos(4)/2+a1t1Offset];
a1t2Pos = [ax1Pos(3)/2-txtSz/1.6+a1t2Offset 10+txtSz/3];
a2t1Pos = [8 ax2Pos(4)/2+a2t1Offset];
a2t2Pos = [ax2Pos(3)/2-txtSz/1.6+a2t2Offset 10+txtSz/3];
a3t1Pos = [8 ax3Pos(4)/2+a3t1Offset];
a3t2Pos = [ax3Pos(3)/2-txtSz/1.6+a3t2Offset 10+txtSz/3];

if ~isempty(a1t1) && handles.volume.imSz1>1 && handles.volume.imSz2>1
    setText(findobj('Tag','IOPt1'),handles.axes.a1,...
        ['''Units'',''points'',''Position'',[' num2str(a1t1Pos) ...
        '],''String'',''' a1t1 ''',''FontSize'',' num2str(txtSz) ','...
        '''FontWeight'',''Bold'',''Color'',[' num2str(txtCol) ...
        '],''EdgeColor'',[' num2str(txtCol) '],''HitTest'',''off'',''Tag'',''IOPt1'''])
else
    delete(findobj('Tag','IOPt1'))
end
if ~isempty(a1t2) && handles.volume.imSz1>1 && handles.volume.imSz2>1
    setText(findobj('Tag','IOPt2'),handles.axes.a1,...
        ['''Units'',''points'',''Position'',[' num2str(a1t2Pos) ...
        '],''String'',''' a1t2 ''',''FontSize'',' num2str(txtSz) ','...
        '''FontWeight'',''Bold'',''Color'',[' num2str(txtCol) ...
        '],''EdgeColor'',[' num2str(txtCol) '],''HitTest'',''off'',''Tag'',''IOPt2'''])
else
    delete(findobj('Tag','IOPt2'))
end
if ~isempty(a2t1) && handles.volume.imSz1>1 && handles.volume.imSz3>1
    setText(findobj('Tag','IOPt3'),handles.axes.a2,...
        ['''Units'',''points'',''Position'',[' num2str(a2t1Pos) ...
        '],''String'',''' a2t1 ''',''FontSize'',' num2str(txtSz) ','...
        '''FontWeight'',''Bold'',''Color'',[' num2str(txtCol) ...
        '],''EdgeColor'',[' num2str(txtCol) '],''HitTest'',''off'',''Tag'',''IOPt3'''])
else
    delete(findobj('Tag','IOPt3'))
end
if ~isempty(a2t2) && handles.volume.imSz1>1 && handles.volume.imSz3>1
    setText(findobj('Tag','IOPt4'),handles.axes.a2,...
        ['''Units'',''points'',''Position'',[' num2str(a2t2Pos) ...
        '],''String'',''' a2t2 ''',''FontSize'',' num2str(txtSz) ','...
        '''FontWeight'',''Bold'',''Color'',[' num2str(txtCol) ...
        '],''EdgeColor'',[' num2str(txtCol) '],''HitTest'',''off'',''Tag'',''IOPt4'''])
else
    delete(findobj('Tag','IOPt4'))
end
if ~isempty(a3t1) && handles.volume.imSz2>1 && handles.volume.imSz3>1
    setText(findobj('Tag','IOPt5'),handles.axes.a3,...
        ['''Units'',''points'',''Position'',[' num2str(a3t1Pos) ...
        '],''String'',''' a3t1 ''',''FontSize'',' num2str(txtSz) ','...
        '''FontWeight'',''Bold'',''Color'',[' num2str(txtCol) ...
        '],''EdgeColor'',[' num2str(txtCol) '],''HitTest'',''off'',''Tag'',''IOPt5'''])
else
    delete(findobj('Tag','IOPt5'))
end
if ~isempty(a3t2) && handles.volume.imSz2>1 && handles.volume.imSz3>1
    setText(findobj('Tag','IOPt6'),handles.axes.a3,...
        ['''Units'',''points'',''Position'',[' num2str(a3t2Pos) ...
        '],''String'',''' a3t2 ''',''FontSize'',' num2str(txtSz) ','...
        '''FontWeight'',''Bold'',''Color'',[' num2str(txtCol) ...
        '],''EdgeColor'',[' num2str(txtCol) '],''HitTest'',''off'',''Tag'',''IOPt6'''])
else
    delete(findobj('Tag','IOPt6'))
end

end

function out = getLabelFromIndex(idx)

switch idx
    case 1
        out = ' F ';
    case 3
        out = ' H ';
    case 5
        out = ' A ';
    case 7
        out = ' P ';
    case 9
        out = ' R ';
    case 11
        out = ' L ';
end

end

function stopDragging(varargin)

global handles

handles.display.crosshairs.crossBoo = true;
toggleCrosshair()
moveArrows()
turnOnArrows()

if ~handles.display.images.fusedImages
    set(handles.figure.f1,'WindowButtonMotionFcn',@figMotionFcn)
    figMotionFcn()
else
    set(handles.figure.f1,'WindowButtonMotionFcn',@figMotionFcnFused)
    figMotionFcnFused()
end

end

function axButtonDownFcnStub(a)

global handles

clickType = get(gcbf,'SelectionType');

if strcmp(clickType,'open')
    if ~handles.display.images.fusedImages && ~strcmpi(handles.volume.modality,'CT')
        invertCM()
    end
    return
end

pt = get(a,'CurrentPoint');
switch(a)
    case handles.axes.a1
        pt_x = round(pt(1,1));
        pt_y = round(pt(1,2));
        % %CONSTRAIN MOUSE CLICKS WITHIN VOLUME RANGE
        % if pt_x < 1 || pt_x > handles.volume.imSz1 || pt_y < 1 || pt_y > handles.volume.imSz2, return; end
    case handles.axes.a2
        pt_x = round(pt(1,1));
        pt_y = handles.volume.imSz3-round(pt(1,2))+1;
        % %CONSTRAIN MOUSE CLICKS WITHIN VOLUME RANGE
        % if pt_x < 1 || pt_x > handles.volume.imSz1 || pt_y < 1 || pt_y > handles.volume.imSz3, return; end
    case handles.axes.a3
        pt_x = round(pt(1,1));
        pt_y = handles.volume.imSz3-round(pt(1,2))+1;
        % %CONSTRAIN MOUSE CLICKS WITHIN VOLUME RANGE
        % if pt_x < 1 || pt_x > handles.volume.imSz2 || pt_y < 1 || pt_y > handles.volume.imSz3, return; end
end

[x,y,z] = getCoordinates(a,pt_x,pt_y);

if handles.display.cursor.overAxes
    switch clickType
        
        case 'normal'
            handles.display.images.x = x;
            handles.display.images.y = y;
            handles.display.images.z = z;
            if ~handles.display.images.fusedImages
                switch(a)
                    case handles.axes.a1
                        showImage3(handles.volume.vol,handles.display.images.h3,handles.display.images.x,handles.display.images.fltr)
                        set(handles.display.sliders.s3,'Value',handles.display.images.x)
                        %SLIDER TEXT
                        axesSliderText(handles.display.sliders.s3t);
                        showImage2(handles.volume.vol,handles.display.images.h2,handles.display.images.y,handles.display.images.fltr)
                        set(handles.display.sliders.s2,'Value',handles.display.images.y)
                        %SLIDER TEXT
                        axesSliderText(handles.display.sliders.s2t);
                        set(handles.figure.f1,'WindowButtonMotionFcn',{@a1DraggingFcn,0})
                    case handles.axes.a2
                        showImage3(handles.volume.vol,handles.display.images.h3,handles.display.images.x,handles.display.images.fltr)
                        set(handles.display.sliders.s3,'Value',handles.display.images.x)
                        %SLIDER TEXT
                        axesSliderText(handles.display.sliders.s3t)
                        showImage1(handles.volume.vol,handles.display.images.h1,handles.display.images.z,handles.display.images.fltr)
                        set(handles.display.sliders.s1,'Value',handles.display.images.z)
                        %SLIDER TEXT
                        axesSliderText(handles.display.sliders.s1t)
                        set(handles.figure.f1,'WindowButtonMotionFcn',{@a2DraggingFcn,0})
                    case handles.axes.a3
                        showImage2(handles.volume.vol,handles.display.images.h2,handles.display.images.y,handles.display.images.fltr)
                        set(handles.display.sliders.s2,'Value',handles.display.images.y)
                        %SLIDER TEXT
                        axesSliderText(handles.display.sliders.s2t);
                        showImage1(handles.volume.vol,handles.display.images.h1,handles.display.images.z,handles.display.images.fltr)
                        set(handles.display.sliders.s1,'Value',handles.display.images.z)
                        %SLIDER TEXT
                        axesSliderText(handles.display.sliders.s1t)
                        set(handles.figure.f1,'WindowButtonMotionFcn',{@a3DraggingFcn,0})
                end
            else
                switch(a)
                    case handles.axes.a1
                        showFusedImage3(handles.axes.a3,handles.volume.vol_1,handles.volume.vol_2,...
                            handles.display.images.h3_1,handles.display.images.h3_2,...
                            handles.display.images.x,handles.display.images.cLim_1,handles.display.images.cLim_2,...
                            handles.display.images.cm_1,handles.display.images.cm_2,handles.display.images.fltr)
                        set(handles.display.sliders.s3,'Value',handles.display.images.x)
                        %SLIDER TEXT
                        axesSliderText(handles.display.sliders.s3t)
                        showFusedImage2(handles.axes.a2,handles.volume.vol_1,handles.volume.vol_2,...
                            handles.display.images.h2_1,handles.display.images.h2_2,...
                            handles.display.images.y,handles.display.images.cLim_1,handles.display.images.cLim_2,...
                            handles.display.images.cm_1,handles.display.images.cm_2,handles.display.images.fltr)
                        set(handles.display.sliders.s2,'Value',handles.display.images.y)
                        %SLIDER TEXT
                        axesSliderText(handles.display.sliders.s2t)
                        set(handles.figure.f1,'WindowButtonMotionFcn',{@a1DraggingFcn,0})
                    case handles.axes.a2
                        showFusedImage3(handles.axes.a3,handles.volume.vol_1,handles.volume.vol_2,...
                            handles.display.images.h3_1,handles.display.images.h3_2,...
                            handles.display.images.x,handles.display.images.cLim_1,handles.display.images.cLim_2,...
                            handles.display.images.cm_1,handles.display.images.cm_2,handles.display.images.fltr)
                        set(handles.display.sliders.s3,'Value',handles.display.images.x)
                        %SLIDER TEXT
                        axesSliderText(handles.display.sliders.s3t)
                        showFusedImage1(handles.axes.a1,handles.volume.vol_1,handles.volume.vol_2,...
                            handles.display.images.h1_1,handles.display.images.h1_2,...
                            handles.display.images.z,handles.display.images.cLim_1,handles.display.images.cLim_2,...
                            handles.display.images.cm_1,handles.display.images.cm_2,handles.display.images.fltr)
                        set(handles.display.sliders.s1,'Value',handles.display.images.z)
                        %SLIDER TEXT
                        axesSliderText(handles.display.sliders.s1t)
                        set(handles.figure.f1,'WindowButtonMotionFcn',{@a2DraggingFcn,0})
                    case handles.axes.a3
                        showFusedImage2(handles.axes.a2,handles.volume.vol_1,handles.volume.vol_2,...
                            handles.display.images.h2_1,handles.display.images.h2_2,...
                            handles.display.images.y,handles.display.images.cLim_1,handles.display.images.cLim_2,...
                            handles.display.images.cm_1,handles.display.images.cm_2,handles.display.images.fltr)
                        set(handles.display.sliders.s2,'Value',handles.display.images.y)
                        %SLIDER TEXT
                        axesSliderText(handles.display.sliders.s2t)
                        showFusedImage1(handles.axes.a1,handles.volume.vol_1,handles.volume.vol_2,...
                            handles.display.images.h1_1,handles.display.images.h1_2,...
                            handles.display.images.z,handles.display.images.cLim_1,handles.display.images.cLim_2,...
                            handles.display.images.cm_1,handles.display.images.cm_2,handles.display.images.fltr)
                        set(handles.display.sliders.s1,'Value',handles.display.images.z)
                        %SLIDER TEXT
                        axesSliderText(handles.display.sliders.s1t)
                        set(handles.figure.f1,'WindowButtonMotionFcn',{@a3DraggingFcn,0})
                end
            end
            turnOnCrosshair()
            %         moveArrows()
            
            if handles.display.images.drawROI
                if handles.mask.draw
                    if handles.mask.mask(x,y,z) %mask point already selected
                        handles.mask.mask(x,y,z) = 0;
                    else
                        handles.mask.mask(x,y,z) = 1;
                    end
                end
                showMask()
                switch(a)
                    case handles.axes.a1
                        set(handles.figure.f1,'WindowButtonMotionFcn',...
                            {@a1DraggingFcn,handles.mask.mask(x,y,z)})
                    case handles.axes.a2
                        set(handles.figure.f1,'WindowButtonMotionFcn',...
                            {@a2DraggingFcn,handles.mask.mask(x,y,z)})
                    case handles.axes.a3
                        set(handles.figure.f1,'WindowButtonMotionFcn',...
                            {@a3DraggingFcn,handles.mask.mask(x,y,z)})
                end
            end
            
        case 'alt'
            %     toggleCrosshair()
            if handles.display.images.drawROI
                if handles.mask.draw
                    switch(a)
                        case handles.axes.a1
                            seedIdx = sub2ind([handles.volume.imSz1,handles.volume.imSz2],x,y);
                            %first check if click was over 1 or 0 pixel in mask
                            CC = bwconncomp(~squeeze(handles.mask.mask(:,:,handles.display.images.z)),4);
                            idx = [];
                            for n=1:CC.NumObjects
                                if ~isempty(find(CC.PixelIdxList{n}==seedIdx,1))
                                    idx = n;
                                end
                            end
                            if ~isempty(idx) %click was over 0 so fill enclosed outline
                                if numel(CC.PixelIdxList)>1 %only fill in connected outline
                                    handles.mask.mask(CC.PixelIdxList{idx}+(z-1)*...
                                        handles.volume.imSz1*handles.volume.imSz2) = 1;
                                end
                            else %click was over 1 so remove connected pixels
                                CC = bwconncomp(squeeze(handles.mask.mask(:,:,handles.display.images.z)),4);
                                for n=1:CC.NumObjects
                                    if ~isempty(find(CC.PixelIdxList{n}==seedIdx,1))
                                        idx = n;
                                    end
                                end
                                handles.mask.mask(CC.PixelIdxList{idx}+(z-1)*...
                                    handles.volume.imSz1*handles.volume.imSz2) = 0;
                            end
                        case handles.axes.a2
                            seedIdx = sub2ind([handles.volume.imSz1,handles.volume.imSz3],x,z);
                            %first check if click was over 1 or 0 pixel in mask
                            CC = bwconncomp(~squeeze(handles.mask.mask(:,handles.display.images.y,:)),4);
                            idx = [];
                            for n=1:CC.NumObjects
                                if ~isempty(find(CC.PixelIdxList{n}==seedIdx,1))
                                    idx = n;
                                end
                            end
                            if ~isempty(idx) %click was over 0 so fill enclosed outline
                                if numel(CC.PixelIdxList)>1 %only fill in connected outline
                                    handles.mask.mask((y-1)*handles.volume.imSz1 + ...
                                        floor((CC.PixelIdxList{idx}-1)/CC.ImageSize(1))*...
                                        handles.volume.imSz1*handles.volume.imSz2 + ...
                                        mod((CC.PixelIdxList{idx}-1),CC.ImageSize(1))+1) = 1;
                                end
                            else %click was over 1 so remove connected pixels
                                CC = bwconncomp(squeeze(handles.mask.mask(:,handles.display.images.y,:)),4);
                                for n=1:CC.NumObjects
                                    if ~isempty(find(CC.PixelIdxList{n}==seedIdx,1))
                                        idx = n;
                                    end
                                end
                                handles.mask.mask((y-1)*handles.volume.imSz1 + ...
                                    floor((CC.PixelIdxList{idx}-1)/CC.ImageSize(1))*...
                                    handles.volume.imSz1*handles.volume.imSz2 + ...
                                    mod((CC.PixelIdxList{idx}-1),CC.ImageSize(1))+1) = 0;
                            end
                        case handles.axes.a3
                            seedIdx = sub2ind([handles.volume.imSz2,handles.volume.imSz3],y,z);
                            %first check if click was over 1 or 0 pixel in mask
                            CC = bwconncomp(~squeeze(handles.mask.mask(handles.display.images.x,:,:)),4);
                            idx = [];
                            for n=1:CC.NumObjects
                                if ~isempty(find(CC.PixelIdxList{n}==seedIdx,1))
                                    idx = n;
                                end
                            end
                            if ~isempty(idx) %click was over 0 so fill enclosed outline
                                if numel(CC.PixelIdxList)>1 %only fill in connected outline
                                    handles.mask.mask(x + ...
                                        floor((CC.PixelIdxList{idx}-1)/CC.ImageSize(1))*...
                                        handles.volume.imSz1*handles.volume.imSz2 + ...
                                        mod(CC.PixelIdxList{idx}-1,CC.ImageSize(1))*handles.volume.imSz1) = 1;
                                end
                            else %click was over 1 so remove connected pixels
                                CC = bwconncomp(squeeze(handles.mask.mask(handles.display.images.x,:,:)),4);
                                for n=1:CC.NumObjects
                                    if ~isempty(find(CC.PixelIdxList{n}==seedIdx,1))
                                        idx = n;
                                    end
                                end
                                handles.mask.mask(x + ...
                                    floor((CC.PixelIdxList{idx}-1)/CC.ImageSize(1))*...
                                    handles.volume.imSz1*handles.volume.imSz2 + ...
                                    mod(CC.PixelIdxList{idx}-1,CC.ImageSize(1))*handles.volume.imSz1) = 0;
                            end
                    end
                    showMask()
                end
            end
    end
end

end

function a1ButtonDownFcn(hObj,~)

axButtonDownFcnStub(hObj)

end

function a2ButtonDownFcn(hObj,~)

axButtonDownFcnStub(hObj)

end

function a3ButtonDownFcn(hObj,~)

axButtonDownFcnStub(hObj)

end

function axDraggingFcnStub(a,drawBool)

global handles

pt = get(a,'currentpoint');

switch a
    case handles.axes.a1
        pt_x = round(pt(1,1));
        pt_y = round(pt(1,2));
    case handles.axes.a2
        pt_x = round(pt(1,1));
        pt_y = handles.volume.imSz3-round(pt(1,2))+1;
    case handles.axes.a3
        pt_x = round(pt(1,1));
        pt_y = handles.volume.imSz3-round(pt(1,2))+1;
end

[x,y,z] = getCoordinates(a,pt_x,pt_y);

handles.display.images.x = x;
handles.display.images.y = y;
handles.display.images.z = z;
if ~handles.display.images.fusedImages
    switch a
        case handles.axes.a1
            showImage3(handles.volume.vol,handles.display.images.h3,handles.display.images.x,handles.display.images.fltr)
            set(handles.display.sliders.s3,'Value',handles.display.images.x)
            %SLIDER TEXT
            axesSliderText(handles.display.sliders.s3t)
            showImage2(handles.volume.vol,handles.display.images.h2,handles.display.images.y,handles.display.images.fltr)
            set(handles.display.sliders.s2,'Value',handles.display.images.y)
            %SLIDER TEXT
            axesSliderText(handles.display.sliders.s2t)
        case handles.axes.a2
            showImage3(handles.volume.vol,handles.display.images.h3,handles.display.images.x,handles.display.images.fltr)
            set(handles.display.sliders.s3,'Value',handles.display.images.x)
            %SLIDER TEXT
            axesSliderText(handles.display.sliders.s3t)
            showImage1(handles.volume.vol,handles.display.images.h1,handles.display.images.z,handles.display.images.fltr)
            set(handles.display.sliders.s1,'Value',handles.display.images.z)
            %SLIDER TEXT
            axesSliderText(handles.display.sliders.s1t)
        case handles.axes.a3
            showImage2(handles.volume.vol,handles.display.images.h2,handles.display.images.y,handles.display.images.fltr)
            set(handles.display.sliders.s2,'Value',handles.display.images.y)
            %SLIDER TEXT
            axesSliderText(handles.display.sliders.s2t)
            showImage1(handles.volume.vol,handles.display.images.h1,handles.display.images.z,handles.display.images.fltr)
            set(handles.display.sliders.s1,'Value',handles.display.images.z)
            %SLIDER TEXT
            axesSliderText(handles.display.sliders.s1t)
    end
    voxelValue(x,y,z)
    coordinatesValue(x,y,z)
else
    switch a
        case handles.axes.a1
            showFusedImage3(handles.axes.a3,handles.volume.vol_1,handles.volume.vol_2,...
                handles.display.images.h3_1,handles.display.images.h3_2,...
                handles.display.images.x,handles.display.images.cLim_1,handles.display.images.cLim_2,...
                handles.display.images.cm_1,handles.display.images.cm_2,handles.display.images.fltr)
            set(handles.display.sliders.s3,'Value',handles.display.images.x)
            %SLIDER TEXT
            axesSliderText(handles.display.sliders.s3t)
            showFusedImage2(handles.axes.a2,handles.volume.vol_1,handles.volume.vol_2,...
                handles.display.images.h2_1,handles.display.images.h2_2,...
                handles.display.images.y,handles.display.images.cLim_1,handles.display.images.cLim_2,...
                handles.display.images.cm_1,handles.display.images.cm_2,handles.display.images.fltr)
            set(handles.display.sliders.s2,'Value',handles.display.images.y)
            %SLIDER TEXT
            axesSliderText(handles.display.sliders.s2t)
        case handles.axes.a2
            showFusedImage3(handles.axes.a3,handles.volume.vol_1,handles.volume.vol_2,...
                handles.display.images.h3_1,handles.display.images.h3_2,...
                handles.display.images.x,handles.display.images.cLim_1,handles.display.images.cLim_2,...
                handles.display.images.cm_1,handles.display.images.cm_2,handles.display.images.fltr)
            set(handles.display.sliders.s3,'Value',handles.display.images.x)
            %SLIDER TEXT
            axesSliderText(handles.display.sliders.s3t)
            showFusedImage1(handles.axes.a1,handles.volume.vol_1,handles.volume.vol_2,...
                handles.display.images.h1_1,handles.display.images.h1_2,...
                handles.display.images.z,handles.display.images.cLim_1,handles.display.images.cLim_2,...
                handles.display.images.cm_1,handles.display.images.cm_2,handles.display.images.fltr)
            set(handles.display.sliders.s1,'Value',handles.display.images.z)
            %SLIDER TEXT
            axesSliderText(handles.display.sliders.s1t)
        case handles.axes.a3
            showFusedImage2(handles.axes.a2,handles.volume.vol_1,handles.volume.vol_2,...
                handles.display.images.h2_1,handles.display.images.h2_2,...
                handles.display.images.y,handles.display.images.cLim_1,handles.display.images.cLim_2,...
                handles.display.images.cm_1,handles.display.images.cm_2,handles.display.images.fltr)
            set(handles.display.sliders.s2,'Value',handles.display.images.y)
            %SLIDER TEXT
            axesSliderText(handles.display.sliders.s2t)
            showFusedImage1(handles.axes.a1,handles.volume.vol_1,handles.volume.vol_2,...
                handles.display.images.h1_1,handles.display.images.h1_2,...
                handles.display.images.z,handles.display.images.cLim_1,handles.display.images.cLim_2,...
                handles.display.images.cm_1,handles.display.images.cm_2,handles.display.images.fltr)
            set(handles.display.sliders.s1,'Value',handles.display.images.z)
            %SLIDER TEXT
            axesSliderText(handles.display.sliders.s1t)
    end
    voxelValueFused(x,y,z)
end

moveCrosshair(a)
% moveArrows()

if handles.display.text.showLocs
    locValue(x,y,z)
end

if handles.display.images.drawROI
    if handles.mask.draw
        if drawBool
            handles.mask.mask(x,y,z) = 1;
        else
            handles.mask.mask(x,y,z) = 0;
        end
    end
    showMask()
end

end

function a1DraggingFcn(~,~,drawBool)

global handles

axDraggingFcnStub(handles.axes.a1,drawBool)

end

function a2DraggingFcn(~,~,drawBool)

global handles

axDraggingFcnStub(handles.axes.a2,drawBool)

end

function a3DraggingFcn(~,~,drawBool)

global handles

axDraggingFcnStub(handles.axes.a3,drawBool)

end

function showMask

global handles

%a1 patch
slice = handles.mask.mask(:,:,handles.display.images.z); %mask indeces
xIdx = handles.mask.a1Xidx(slice)';
yIdx = handles.mask.a1Yidx(slice)';
xD = [xIdx-0.5;xIdx+0.5;xIdx+0.5;xIdx-0.5;xIdx-0.5];
yD = handles.volume.imSz2+1-[yIdx-0.5;yIdx-0.5;yIdx+0.5;yIdx+0.5;yIdx-0.5];
set(handles.mask.a1patch,'Xdata',xD,'Ydata',yD)
%a2 patch
slice = handles.mask.mask(:,handles.display.images.y,:); %mask indeces
xIdx = handles.mask.a2Xidx(slice)';
yIdx = handles.mask.a2Yidx(slice)';
xD = [xIdx-0.5;xIdx+0.5;xIdx+0.5;xIdx-0.5;xIdx-0.5];
yD = [yIdx-0.5;yIdx-0.5;yIdx+0.5;yIdx+0.5;yIdx-0.5];
set(handles.mask.a2patch,'Xdata',xD,'Ydata',yD)
%a3 patch
slice = handles.mask.mask(handles.display.images.x,:,:); %mask indeces
xIdx = handles.mask.a3Xidx(slice)';
yIdx = handles.mask.a3Yidx(slice)';
xD = [xIdx-0.5;xIdx+0.5;xIdx+0.5;xIdx-0.5;xIdx-0.5];
yD = [yIdx-0.5;yIdx-0.5;yIdx+0.5;yIdx+0.5;yIdx-0.5];
set(handles.mask.a3patch,'Xdata',xD,'Ydata',yD)

end

function defineDrawButton

global handles

figPos = get(handles.figure.f1,'Position');
figSz = figPos(3:4);

decBtnWidth = 90/handles.figure.DPIscl;
decBtnHeight = 36/handles.figure.DPIscl;
decBtn_y = 0.05*figSz(2);
decBtn1_x = 0.8*figPos(3);

doneBtnPos = [decBtn1_x-decBtnWidth/2 decBtn_y-decBtnHeight/2 decBtnWidth decBtnHeight];

fs = 0.006*figSz(1)*handles.figure.DPIscl;

str = ['''Position'',[' num2str(doneBtnPos) '],'...
    '''FontSize'',' num2str(fs) ','...
    '''HorizontalAlignment'',''Center'','...
    '''Enable'',''Off'','...
    '''String'',''Draw ROI'','...
    '''BackgroundColor'',[' num2str(get(findobj('Tag','3DviewerFig'),'color')) ']'];
handles.drawBtn = initiateObject('drawBtn','uicontrol',str);

end

function defineDoneButton

global handles

figPos = get(handles.figure.f1,'Position');
figSz = figPos(3:4);

decBtnWidth = 90/handles.figure.DPIscl;
decBtnHeight = 36/handles.figure.DPIscl;
decBtn_y = 0.05*figSz(2);
decBtn1_x = 0.8*figPos(3) + (decBtn_y-decBtnHeight/2) + 0.05*figSz(1);

doneBtnPos = [decBtn1_x-decBtnWidth/2 decBtn_y-decBtnHeight/2 decBtnWidth decBtnHeight];

fs = 0.006*figSz(1)*handles.figure.DPIscl;

str = ['''Position'',[' num2str(doneBtnPos) '],'...
    '''FontSize'',' num2str(fs) ','...
    '''HorizontalAlignment'',''Center'','...
    '''Enable'',''Off'','...
    '''String'',''Done'','...
    '''BackgroundColor'',[' num2str(get(findobj('Tag','3DviewerFig'),'color')) ']'];
handles.doneBtn = initiateObject('doneBtn','uicontrol',str);

end

function drawButtonFcn(hObj,~)

global handles

if handles.mask.draw
    handles.mask.draw = false;
    if handles.display.text.instructionsJavaBool
        handles.display.text.instructions.Text = [handles.display.text.instructionsStr{1:end-1} '</html>'];
    else
        set(handles.display.text.instructions,'String',handles.display.text.instructionsStr(1:end-1))
    end
    set(hObj,'BackgroundColor',get(handles.figure.f1,'Color'))
else
    handles.mask.draw = true;
    if handles.display.text.instructionsJavaBool
        handles.display.text.instructions.Text = [handles.display.text.instructionsStr{:} '</html>'];
    else
        set(handles.display.text.instructions,'String',handles.display.text.instructionsStr)
    end
    set(hObj,'BackgroundColor',[0.4660 0.6740 0.1880])
end

end

function setText(h,a,str)

if ~isempty(h)
    eval(['set(h,' str ');'])
else
    axes(a)
    eval(['text(' str ');'])
end

end

function invertCM(varargin)

setCM(1-colormap)

end

function setCM(cm)

colormap(cm)
% set(findobj('Type','axes'),'Color',cm(1,:))
set(findobj('Tag','axes1','-or','Tag','axes2','-or','Tag','axes3'),'Color',cm(1,:))

end

function CMselect(h,eventdata)

global handles

if ~handles.display.images.fusedImages
    switch lower(get(eventdata.NewValue,'String'))
        case 'gray'
            if strcmpi(handles.volume.modality,'PET')
                setCM(1-gray)
            else
                setCM(gray)
            end
        case 'hot'
            setCM(hotMetal)
        case 'pet'
            setCM(PET)
    end
else
    switch get(h,'Tag')
        case 'CMPanel_1'
            switch lower(get(eventdata.NewValue,'String'))
                case 'gray'
                    if strcmpi(handles.volume.modality_1,'PET')
                        handles.display.images.cm_1 = 1-gray;
                    else
                        handles.display.images.cm_1 = gray;
                    end
                case 'hot'
                    handles.display.images.cm_1 = hotMetal;
                case 'pet'
                    handles.display.images.cm_1 = PET;
            end
        case 'CMPanel_2'
            switch lower(get(eventdata.NewValue,'String'))
                case 'gray'
                    if strcmpi(handles.volume.modality_2,'PET')
                        handles.display.images.cm_2 = 1-gray;
                    else
                        handles.display.images.cm_2 = gray;
                    end
                case 'hot'
                    handles.display.images.cm_2 = hotMetal;
                case 'pet'
                    handles.display.images.cm_2 = PET;
            end
    end
    setFusedDisplayWindow()
end

end

function CTselect(h,eventdata)

global handles

if ~handles.display.images.fusedImages
    handles.display.images.cLim = getCTwindow(get(eventdata.NewValue,'String'));
    if handles.display.images.cLim(1)<get(findobj('Tag','s5Slider'),'Min')
        handles.display.images.cLim(1) = get(findobj('Tag','s5Slider'),'Min');
    end
    if handles.display.images.cLim(2)>get(findobj('Tag','s4Slider'),'Max')
        handles.display.images.cLim(2) = get(findobj('Tag','s4Slider'),'Max');
    end
    
    set(findobj('Tag','s4Slider'),'Value',handles.display.images.cLim(2))
    set(findobj('Tag','s5Slider'),'Value',handles.display.images.cLim(1))
    
    s4tString = formatValueString(handles.display.images.cLim(2)/handles.display.images.sclFctr);
    s5tString = formatValueString(handles.display.images.cLim(1)/handles.display.images.sclFctr);
    
    set(handles.display.sliders.s4t,'String',s4tString)
    set(handles.display.sliders.s5t,'String',s5tString)
    
    setDisplayWindow()
else
    switch get(h,'Tag')
        case 'CMPanel_1'
            handles.display.images.cLim_1 = getCTwindow(get(eventdata.NewValue,'String'));
        case 'CMPanel_2'
            handles.display.images.cLim_2 = getCTwindow(get(eventdata.NewValue,'String'));
    end
    if handles.display.images.cLim_1(1)<get(findobj('Tag','s5Slider_1'),'Min')
        handles.display.images.cLim_1(1) = get(findobj('Tag','s5Slider_1'),'Min');
    end
    if handles.display.images.cLim_1(2)>get(findobj('Tag','s4Slider_1'),'Max')
        handles.display.images.cLim_1(2) = get(findobj('Tag','s4Slider_1'),'Max');
    end
    if handles.display.images.cLim_2(1)<get(findobj('Tag','s5Slider_2'),'Min')
        handles.display.images.cLim_2(1) = get(findobj('Tag','s5Slider_2'),'Min');
    end
    if handles.display.images.cLim_2(2)>get(findobj('Tag','s4Slider_2'),'Max')
        handles.display.images.cLim_2(2) = get(findobj('Tag','s4Slider_2'),'Max');
    end
    set(findobj('Tag','s4Slider_1'),'Value',handles.display.images.cLim_1(2))
    set(findobj('Tag','s5Slider_1'),'Value',handles.display.images.cLim_1(1))
    set(findobj('Tag','s4Slider_2'),'Value',handles.display.images.cLim_2(2))
    set(findobj('Tag','s5Slider_2'),'Value',handles.display.images.cLim_2(1))
    
    s4tString_1 = formatValueString(handles.display.images.cLim_1(2)/handles.display.images.sclFctr_1);
    s5tString_1 = formatValueString(handles.display.images.cLim_1(1)/handles.display.images.sclFctr_1);
    s4tString_2 = formatValueString(handles.display.images.cLim_2(2)/handles.display.images.sclFctr_2);
    s5tString_2 = formatValueString(handles.display.images.cLim_2(1)/handles.display.images.sclFctr_2);
    
    set(handles.display.sliders.s4t_1,'String',s4tString_1)
    set(handles.display.sliders.s5t_1,'String',s5tString_1)
    set(handles.display.sliders.s4t_2,'String',s4tString_2)
    set(handles.display.sliders.s5t_2,'String',s5tString_2)
    
    setFusedDisplayWindow()
end

end

function cLim = getCTwindow(windowName)

switch lower(windowName)
    case 'body'
        cLim = [40-200 40+200];
    case 'lung'
        cLim = [-600-600 -600+600];
    case 'bone'
        cLim = [450-750 450+750];
    otherwise
        error(['CT window not defined for ' windowName])
end

end

function hw = getSliderKnobWidth(sWidth,arrowButtonWidth,ss)

if numel(ss)==2, ss = ss(2); end

sWidth = sWidth-2*arrowButtonWidth;

pC = [0.0212 -0.2772 0.7414 0.0116]; %FOUND EMPIRICALLY

hw = sWidth * (pC(1)*ss^3+pC(2)*ss^2+pC(3)*ss+pC(4));

end

function axesSliderText(h)

global handles

switch h
    case (handles.display.sliders.s1t)
        value = handles.display.images.z;
        pos = [handles.display.sliders.s1tLocs(handles.display.images.z-(handles.display.images.zLimits(1)-1)) ...
            handles.display.sliders.st_y handles.display.sliders.stSz];
    case (handles.display.sliders.s2t)
        value = handles.display.images.y;
        pos = [handles.display.sliders.s2tLocs(handles.display.images.y-(handles.display.images.yLimits(1)-1)) ...
            handles.display.sliders.st_y handles.display.sliders.stSz];
    case (handles.display.sliders.s3t)
        value = handles.display.images.x;
        pos = [handles.display.sliders.s3tLocs(handles.display.images.x-(handles.display.images.xLimits(1)-1)) ...
            handles.display.sliders.st_y handles.display.sliders.stSz];
end

% set(h,'Position',pos,'String',num2str(value))
set(h,'String',num2str(value))

end

function zoomFcn(~,eventdata)

global handles

maxZoomFctr = 6;

if ~(handles.display.images.zoomFctr==maxZoomFctr && eventdata.VerticalScrollCount<0) || ...
        ~(handles.display.images.zoomFctr==1 && eventdata.VerticalScrollCount>0)
    
    handles.display.images.zoomFctr = handles.display.images.zoomFctr-eventdata.VerticalScrollCount;
    if handles.display.images.zoomFctr<1, handles.display.images.zoomFctr = 1; end
    if handles.display.images.zoomFctr>maxZoomFctr, handles.display.images.zoomFctr = maxZoomFctr; end
    
    updateZoomWindow()
end

end

function updateZoomWindow

global handles

currentrange1 = (handles.volume.imSz1-1)/handles.display.images.zoomFctr;
currentrange2 = (handles.volume.imSz2-1)/handles.display.images.zoomFctr;
currentrange3 = (handles.volume.imSz3-1)/handles.display.images.zoomFctr;
handles.display.images.xLimits = round([handles.display.images.x-currentrange1/2 handles.display.images.x+currentrange1/2]);
handles.display.images.yLimits = round([handles.display.images.y-currentrange2/2 handles.display.images.y+currentrange2/2]);
handles.display.images.zLimits = round([handles.display.images.z-currentrange3/2 handles.display.images.z+currentrange3/2]);

%MUST ACCOUNT FOR FLIPPED Y DIMENSION IN FIRST AXIS
a1yLimits = round([handles.volume.imSz2-handles.display.images.y+1-currentrange2/2 ...
    handles.volume.imSz2-handles.display.images.y+1+currentrange2/2]);

if handles.display.images.xLimits(1)<1
    handles.display.images.xLimits = [1 round(1+currentrange1)];
elseif handles.display.images.xLimits(2)>handles.volume.imSz1
    handles.display.images.xLimits = [round(handles.volume.imSz1-currentrange1) handles.volume.imSz1];
end
if handles.display.images.yLimits(1)<1
    handles.display.images.yLimits = [1 round(1+currentrange2)];
elseif handles.display.images.yLimits(2)>handles.volume.imSz2
    handles.display.images.yLimits = [round(handles.volume.imSz2-currentrange2) handles.volume.imSz2];
end
if handles.display.images.zLimits(1)<1
    handles.display.images.zLimits = [1 round(1+currentrange3)];
elseif handles.display.images.zLimits(2)>handles.volume.imSz3
    handles.display.images.zLimits = [round(handles.volume.imSz3-currentrange3) handles.volume.imSz3];
end
if a1yLimits(1)<1
    a1yLimits = [1 round(1+currentrange2)];
elseif a1yLimits(2)>handles.volume.imSz2
    a1yLimits = [round(handles.volume.imSz2-currentrange2) handles.volume.imSz2];
end

a1XLim = handles.display.images.xLimits+handles.display.images.xlimPad/handles.display.images.zoomFctr;
a1YLim = a1yLimits+handles.display.images.ylimPad/handles.display.images.zoomFctr;
a2XLim = handles.display.images.xLimits+handles.display.images.xlimPad/handles.display.images.zoomFctr;
a2YLim = handles.display.images.zLimits+handles.display.images.ylimPad/handles.display.images.zoomFctr;
a3XLim = handles.display.images.yLimits+handles.display.images.xlimPad/handles.display.images.zoomFctr;
a3YLim = handles.display.images.zLimits+handles.display.images.ylimPad/handles.display.images.zoomFctr;

set(handles.axes.a1,'XLim',a1XLim,'Ylim',a1YLim)
set(handles.axes.a2,'XLim',a2XLim,'Ylim',a2YLim)
set(handles.axes.a3,'XLim',a3XLim,'Ylim',a3YLim)

%REDEFINE ARROW LOCATIONS
initializeArrows()

if ~handles.display.images.fusedImages
    figMotionFcn()
else
    figMotionFcnFused()
end

ax1Pos = get(handles.axes.a1,'Position');
ax2Pos = get(handles.axes.a2,'Position');
ax3Pos = get(handles.axes.a3,'Position');
s1Pos = get(handles.display.sliders.s1,'Position');
s2Pos = get(handles.display.sliders.s2,'Position');
s3Pos = get(handles.display.sliders.s3,'Position');
arrowButtonWidth = 16; %TO ACCOUNT FOR SLIDER "TROUGH" SIZE

%REDEFINE NAVIGATION SLIDERS AND COMPUTE ARRRAYS OF TEXT BOX LOCATIONS
if handles.volume.imSz3>1
    set(handles.display.sliders.s1,'Min',handles.display.images.zLimits(1),'Max',handles.display.images.zLimits(2),...
        'SliderStep',[1/diff(handles.display.images.zLimits) 5/diff(handles.display.images.zLimits)])
    shw = getSliderKnobWidth(ax1Pos(3),arrowButtonWidth,get(handles.display.sliders.s1,'SliderStep'));
    sm = (s1Pos(3)-shw-2*arrowButtonWidth)/diff(handles.display.images.zLimits);
    sb = s1Pos(1)+arrowButtonWidth+shw/2-...
        ((s1Pos(3)-2*arrowButtonWidth-shw)/diff(handles.display.images.zLimits));
    handles.display.sliders.s1tLocs = ...
        sm.*(1:diff(handles.display.images.zLimits)+1)+sb-handles.display.sliders.stSz(1)/2;
else
    set(handles.display.sliders.s1,'Min',handles.display.images.zLimits(1),'Max',handles.display.images.zLimits(2))
    handles.display.sliders.s1tLocs = s1Pos(1)+s1Pos(3)/2-handles.display.sliders.stSz(1)/2;
end
if handles.volume.imSz2>1
    set(handles.display.sliders.s2,'Min',handles.display.images.yLimits(1),'Max',handles.display.images.yLimits(2),...
        'SliderStep',[1/diff(handles.display.images.yLimits) 5/diff(handles.display.images.yLimits)])
    shw = getSliderKnobWidth(ax2Pos(3),arrowButtonWidth,get(handles.display.sliders.s2,'SliderStep'));
    sm = (s2Pos(3)-shw-2*arrowButtonWidth)/diff(handles.display.images.yLimits);
    sb = s2Pos(1)+arrowButtonWidth+shw/2-...
        ((s2Pos(3)-2*arrowButtonWidth-shw)/diff(handles.display.images.yLimits));
    handles.display.sliders.s2tLocs = ...
        sm.*(1:diff(handles.display.images.yLimits)+1)+sb-handles.display.sliders.stSz(1)/2;
else
    set(handles.display.sliders.s2,'Min',handles.display.images.yLimits(1),'Max',handles.display.images.yLimits(2))
    handles.display.sliders.s2tLocs = s2Pos(1)+s2Pos(3)/2-handles.display.sliders.stSz(1)/2;
end
if handles.volume.imSz1>1
    set(handles.display.sliders.s3,'Min',handles.display.images.xLimits(1),'Max',handles.display.images.xLimits(2),...
        'SliderStep',[1/diff(handles.display.images.xLimits) 5/diff(handles.display.images.xLimits)])
    shw = getSliderKnobWidth(ax3Pos(3),arrowButtonWidth,get(handles.display.sliders.s3,'SliderStep'));
    sm = (s3Pos(3)-shw-2*arrowButtonWidth)/diff(handles.display.images.xLimits);
    sb = s3Pos(1)+arrowButtonWidth+shw/2-...
        ((s3Pos(3)-2*arrowButtonWidth-shw)/diff(handles.display.images.xLimits));
    handles.display.sliders.s3tLocs = ...
        sm.*(1:diff(handles.display.images.xLimits)+1)+sb-handles.display.sliders.stSz(1)/2;
else
    set(handles.display.sliders.s3,'Min',handles.display.images.xLimits(1),'Max',handles.display.images.xLimits(2))
    handles.display.sliders.s2tLocs = s2Pos(1)+s2Pos(3)/2-handles.display.sliders.stSz(1)/2;
end

%ADJUST TEXT BOX POSITIONS (NOT NEEDED IF STATIONARY)
axesSliderText(handles.display.sliders.s1t);
axesSliderText(handles.display.sliders.s2t);
axesSliderText(handles.display.sliders.s3t);

end

function stampImageInfo

global handles

figSz = get(handles.figure.f1,'Position');
figSz = [figSz(3) figSz(4)];

if ~strcmpi(handles.volume.filePath,'{unknown}')
    txtStr1 = '[{''image file path: ''}; {''volume matrix size: ''}; {''voxel dimensions (mm): ''}]';
    txtStr2 = ['{''' handles.volume.filePath '''}'];
else
    txtStr1 = '[{''volume matrix size: ''}; {''voxel dimensions (mm): ''}]';
    txtStr2 = '';
end
txtStr2 = [txtStr2 ...
    '; {''' num2str(handles.volume.imSz1) ' x ' num2str(handles.volume.imSz2) ' x ' num2str(handles.volume.imSz3) '''}'];
txtStr2 = [txtStr2 ...
    '; {''' num2str(handles.volume.voxDimX) ' x ' num2str(handles.volume.voxDimY) ' x ' num2str(handles.volume.voxDimZ) '''}'];
txtStr2 = ['[' txtStr2 ']'];


t1Sz = [0.09*figSz(1) 0.05*figSz(2)];
t1Pos = [0.01*figSz(1) figSz(2)-0.011*figSz(2)-t1Sz(2) t1Sz];
t2Sz = [0.8*figSz(1) 0.05*figSz(2)];
t2Pos = [t1Pos(1)+t1Sz(1) figSz(2)-0.011*figSz(2)-t2Sz(2) t2Sz];

str = ['''Style'',''Text'','...
    '''Units'',''Pixel'','...
    '''Position'',[' num2str(t1Pos) '],'...
    '''FontSize'',' num2str(0.005*figSz(1)*handles.figure.DPIscl) ','...
    '''HorizontalAlignment'',''Right'','...
    '''FontWeight'',''Normal'','...
    '''String'',' txtStr1 ','...
    '''BackgroundColor'',[' num2str(get(findobj('Tag','3DviewerFig'),'color')) ']'];
h1 = initiateObject('StampText1','uicontrol',str);

str = ['''Style'',''Text'','...
    '''Units'',''Pixel'','...
    '''Position'',[' num2str(t2Pos) '],'...
    '''FontSize'',' num2str(0.005*figSz(1)*handles.figure.DPIscl) ','...
    '''HorizontalAlignment'',''Left'','...
    '''FontWeight'',''Normal'','...
    '''String'',' txtStr2 ','...
    '''BackgroundColor'',[' num2str(get(findobj('Tag','3DviewerFig'),'color')) ']'];
h2 = initiateObject('StampText2','uicontrol',str);

uistack(h1,'bottom')
uistack(h2,'bottom')

end

function freezeImageColors(a,h)

%converts object color indexes into colormap to true-color data

%adapted from
%  'freezeColors / unfreezeColors' written by John Iversen

%current colormap
cmap = colormap;
nColors = size(cmap,1);
cax = caxis(a);
%  current colormap
g = get(h);
%preserve parent axis clim
originalClim = get(a, 'clim');
%get image CData
cdata = g.CData;
%save original indexed data for use with unfreezeImageColors (not used)
siz = size(cdata);
%convert cdata to indexes into colormap
idx = ceil( (double(cdata) - cax(1)) / (cax(2)-cax(1)) * nColors);
%clamp to [1, nColors]
idx(idx<1) = 1;
idx(idx>nColors) = nColors;
%make true-color data--using current colormap
realcolor = zeros([siz 3]);
for i = 1:3,
    c = cmap(idx,i);
    c = reshape(c,siz);
    realcolor(:,:,i) = c;
end
%replace original CData with true-color data
set(h,'CData',realcolor);
%restore clim (so colorbar will show correct limits)
set(a,'clim',originalClim)

end

function outStr = formatValueString(val)

absval = abs(val);

if absval>=1e4
    %very high
    outStr = sprintf('%0.1e',val);
elseif val==round(val)
    %integer
    if val==0
        outStr = sprintf('%0.0f',0);
    else
        outStr = sprintf('%0.0f',val);
    end
elseif absval<1e-2
    %very low value
    outStr = sprintf('%0.1e',val);
elseif absval<1
    outStr = sprintf('%0.3f',val);
else
    if mod(round(val*1e2)/1e2,1)==0
        outStr = sprintf('%0.0f',val);
    else
        outStr = sprintf('%0.2f',val);
    end
end

end

function h = initiateObject(tag, obj, str, varargin)

%varargin is used to pass button group handle if needed

eval(['test = findobj(''Tag'','''  tag ''');'])
if ~isempty(test)
    eval(['set(test,' str ')'])
    h = test;
else
    if isempty(varargin)
        eval(['h = '  obj '(' str ',''Tag'',''' tag ''');'])
    else
        eval(['h = '  obj '(varargin{1},' str ',''Tag'',''' tag ''');'])
    end
end

end

function cm = hotMetal

%according to http://medical.nema.org/dicom/2013/output/chtml/part06/chapter_B.html

cm = [
    0	0	0
    2	0	0
    4	0	0
    6	0	0
    8	0	0
    10	0	0
    12	0	0
    14	0	0
    16	0	0
    18	0	0
    20	0	0
    22	0	0
    24	0	0
    26	0	0
    28	0	0
    30	0	0
    32	0	0
    34	0	0
    36	0	0
    38	0	0
    40	0	0
    42	0	0
    44	0	0
    46	0	0
    48	0	0
    50	0	0
    52	0	0
    54	0	0
    56	0	0
    58	0	0
    60	0	0
    62	0	0
    64	0	0
    66	0	0
    68	0	0
    70	0	0
    72	0	0
    74	0	0
    76	0	0
    78	0	0
    80	0	0
    82	0	0
    84	0	0
    86	0	0
    88	0	0
    90	0	0
    92	0	0
    94	0	0
    96	0	0
    98	0	0
    100	0	0
    102	0	0
    104	0	0
    106	0	0
    108	0	0
    110	0	0
    112	0	0
    114	0	0
    116	0	0
    118	0	0
    120	0	0
    122	0	0
    124	0	0
    126	0	0
    128	0	0
    130	0	0
    132	0	0
    134	0	0
    136	0	0
    138	0	0
    140	0	0
    142	0	0
    144	0	0
    146	0	0
    148	0	0
    150	0	0
    152	0	0
    154	0	0
    156	0	0
    158	0	0
    160	0	0
    162	0	0
    164	0	0
    166	0	0
    168	0	0
    170	0	0
    172	0	0
    174	0	0
    176	0	0
    178	0	0
    180	0	0
    182	0	0
    184	0	0
    186	0	0
    188	0	0
    190	0	0
    192	0	0
    194	0	0
    196	0	0
    198	0	0
    200	0	0
    202	0	0
    204	0	0
    206	0	0
    208	0	0
    210	0	0
    212	0	0
    214	0	0
    216	0	0
    218	0	0
    220	0	0
    222	0	0
    224	0	0
    226	0	0
    228	0	0
    230	0	0
    232	0	0
    234	0	0
    236	0	0
    238	0	0
    240	0	0
    242	0	0
    244	0	0
    246	0	0
    248	0	0
    250	0	0
    252	0	0
    254	0	0
    255	0	0
    255	2	0
    255	4	0
    255	6	0
    255	8	0
    255	10	0
    255	12	0
    255	14	0
    255	16	0
    255	18	0
    255	20	0
    255	22	0
    255	24	0
    255	26	0
    255	28	0
    255	30	0
    255	32	0
    255	34	0
    255	36	0
    255	38	0
    255	40	0
    255	42	0
    255	44	0
    255	46	0
    255	48	0
    255	50	0
    255	52	0
    255	54	0
    255	56	0
    255	58	0
    255	60	0
    255	62	0
    255	64	0
    255	66	0
    255	68	0
    255	70	0
    255	72	0
    255	74	0
    255	76	0
    255	78	0
    255	80	0
    255	82	0
    255	84	0
    255	86	0
    255	88	0
    255	90	0
    255	92	0
    255	94	0
    255	96	0
    255	98	0
    255	100	0
    255	102	0
    255	104	0
    255	106	0
    255	108	0
    255	110	0
    255	112	0
    255	114	0
    255	116	0
    255	118	0
    255	120	0
    255	122	0
    255	124	0
    255	126	0
    255	128	4
    255	130	8
    255	132	12
    255	134	16
    255	136	20
    255	138	24
    255	140	28
    255	142	32
    255	144	36
    255	146	40
    255	148	44
    255	150	48
    255	152	52
    255	154	56
    255	156	60
    255	158	64
    255	160	68
    255	162	72
    255	164	76
    255	166	80
    255	168	84
    255	170	88
    255	172	92
    255	174	96
    255	176	100
    255	178	104
    255	180	108
    255	182	112
    255	184	116
    255	186	120
    255	188	124
    255	190	128
    255	192	132
    255	194	136
    255	196	140
    255	198	144
    255	200	148
    255	202	152
    255	204	156
    255	206	160
    255	208	164
    255	210	168
    255	212	172
    255	214	176
    255	216	180
    255	218	184
    255	220	188
    255	222	192
    255	224	196
    255	226	200
    255	228	204
    255	230	208
    255	232	212
    255	234	216
    255	236	220
    255	238	224
    255	240	228
    255	242	232
    255	244	236
    255	246	240
    255	248	244
    255	250	248
    255	252	252
    255	255	255
    ] / 255;

end

function cm = PET

%according to http://medical.nema.org/dicom/2013/output/chtml/part06/chapter_B.html

cm = [
    0	0	0
    0	2	1
    0	4	3
    0	6	5
    0	8	7
    0	10	9
    0	12	11
    0	14	13
    0	16	15
    0	18	17
    0	20	19
    0	22	21
    0	24	23
    0	26	25
    0	28	27
    0	30	29
    0	32	31
    0	34	33
    0	36	35
    0	38	37
    0	40	39
    0	42	41
    0	44	43
    0	46	45
    0	48	47
    0	50	49
    0	52	51
    0	54	53
    0	56	55
    0	58	57
    0	60	59
    0	62	61
    0	65	63
    0	67	65
    0	69	67
    0	71	69
    0	73	71
    0	75	73
    0	77	75
    0	79	77
    0	81	79
    0	83	81
    0	85	83
    0	87	85
    0	89	87
    0	91	89
    0	93	91
    0	95	93
    0	97	95
    0	99	97
    0	101	99
    0	103	101
    0	105	103
    0	107	105
    0	109	107
    0	111	109
    0	113	111
    0	115	113
    0	117	115
    0	119	117
    0	121	119
    0	123	121
    0	125	123
    0	128	125
    1	126	127
    3	124	129
    5	122	131
    7	120	133
    9	118	135
    11	116	137
    13	114	139
    15	112	141
    17	110	143
    19	108	145
    21	106	147
    23	104	149
    25	102	151
    27	100	153
    29	98	155
    31	96	157
    33	94	159
    35	92	161
    37	90	163
    39	88	165
    41	86	167
    43	84	169
    45	82	171
    47	80	173
    49	78	175
    51	76	177
    53	74	179
    55	72	181
    57	70	183
    59	68	185
    61	66	187
    63	64	189
    65	63	191
    67	61	193
    69	59	195
    71	57	197
    73	55	199
    75	53	201
    77	51	203
    79	49	205
    81	47	207
    83	45	209
    85	43	211
    86	41	213
    88	39	215
    90	37	217
    92	35	219
    94	33	221
    96	31	223
    98	29	225
    100	27	227
    102	25	229
    104	23	231
    106	21	233
    108	19	235
    110	17	237
    112	15	239
    114	13	241
    116	11	243
    118	9	245
    120	7	247
    122	5	249
    124	3	251
    126	1	253
    128	0	255
    130	2	252
    132	4	248
    134	6	244
    136	8	240
    138	10	236
    140	12	232
    142	14	228
    144	16	224
    146	18	220
    148	20	216
    150	22	212
    152	24	208
    154	26	204
    156	28	200
    158	30	196
    160	32	192
    162	34	188
    164	36	184
    166	38	180
    168	40	176
    170	42	172
    171	44	168
    173	46	164
    175	48	160
    177	50	156
    179	52	152
    181	54	148
    183	56	144
    185	58	140
    187	60	136
    189	62	132
    191	64	128
    193	66	124
    195	68	120
    197	70	116
    199	72	112
    201	74	108
    203	76	104
    205	78	100
    207	80	96
    209	82	92
    211	84	88
    213	86	84
    215	88	80
    217	90	76
    219	92	72
    221	94	68
    223	96	64
    225	98	60
    227	100	56
    229	102	52
    231	104	48
    233	106	44
    235	108	40
    237	110	36
    239	112	32
    241	114	28
    243	116	24
    245	118	20
    247	120	16
    249	122	12
    251	124	8
    253	126	4
    255	128	0
    255	130	4
    255	132	8
    255	134	12
    255	136	16
    255	138	20
    255	140	24
    255	142	28
    255	144	32
    255	146	36
    255	148	40
    255	150	44
    255	152	48
    255	154	52
    255	156	56
    255	158	60
    255	160	64
    255	162	68
    255	164	72
    255	166	76
    255	168	80
    255	170	85
    255	172	89
    255	174	93
    255	176	97
    255	178	101
    255	180	105
    255	182	109
    255	184	113
    255	186	117
    255	188	121
    255	190	125
    255	192	129
    255	194	133
    255	196	137
    255	198	141
    255	200	145
    255	202	149
    255	204	153
    255	206	157
    255	208	161
    255	210	165
    255	212	170
    255	214	174
    255	216	178
    255	218	182
    255	220	186
    255	222	190
    255	224	194
    255	226	198
    255	228	202
    255	230	206
    255	232	210
    255	234	214
    255	236	218
    255	238	222
    255	240	226
    255	242	230
    255	244	234
    255	246	238
    255	248	242
    255	250	246
    255	252	250
    255	255	255
    ] / 255;

end