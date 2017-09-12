function [bedRangeX, bedRangeY, bedRangeZ] = getScannerBedRanges(varargin)

%this function gets the bed coverage from IOP, IPP, image matrix size, and pixel spacing
%the output is in SCANNER coordinates, not image coordinates unless IOP==[1 0 0 0 1 0]

if nargin==1 || nargin==7
    [IOP, IPP, imSz1, imSz2, voxDim1, voxDim2, voxDim3] = parseInputs(varargin{:});
elseif nargin==6
    [IOP, IPP, imSz1, imSz2, voxDim1, voxDim2] = parseInputs(varargin{:});
end
imSz3 = size(IPP,2);
if ~exist('voxDim3','var')
    if size(IPP,2)>1
        voxDim3 = norm(IPP(:,2)-IPP(:,1));
    else %single slice
        error('Must input thickness of a single slice')
    end
end

voxDim1 = abs(voxDim1);
voxDim2 = abs(voxDim2);
voxDim3 = abs(voxDim3);

%calculate all voxel positions
[locX, locY, locZ] = getVoxelLocs(IOP, IPP, imSz1, imSz2, imSz3, voxDim1, voxDim2);
%find voxel dimensions along scanner coordinate
[voxDimX, voxDimY, voxDimZ] = rotateValues(voxDim1,voxDim2,voxDim3,IOP);

if locX(1,1,1)<locX(end,end,end)
    bedPtsX = [min(locX(:)) max(locX(:))];
elseif locX(1,1,1)>locX(end,end,end)
    bedPtsX = [max(locX(:)) min(locX(:))];
else
    error('Problem here')
end
if locY(1,1,1)<locY(end,end,end)
    bedPtsY = [min(locY(:)) max(locY(:))];
elseif locY(1,1,1)>locY(end,end,end)
    bedPtsY = [max(locY(:)) min(locY(:))];
else
    error('Problem here')
end
if locZ(1,1,1)<locZ(end,end,end)
    bedPtsZ = [min(locZ(:)) max(locZ(:))];
elseif locZ(1,1,1)>locZ(end,end,end)
    bedPtsZ = [max(locZ(:)) min(locZ(:))];
else
    error('Problem here')
end
%define bed ranges
if bedPtsX(1)<=bedPtsX(2)
    bedRangeX = [bedPtsX(1)-abs(voxDimX)/2 bedPtsX(2)+abs(voxDimX)/2];
else
    bedRangeX = [bedPtsX(1)+abs(voxDimX)/2 bedPtsX(2)-abs(voxDimX)/2];
end
if bedPtsY(1)<=bedPtsY(2)
    bedRangeY = [bedPtsY(1)-abs(voxDimY)/2 bedPtsY(2)+abs(voxDimY)/2];
else
    bedRangeY = [bedPtsY(1)+abs(voxDimY)/2 bedPtsY(2)-abs(voxDimY)/2];
end
if bedPtsZ(1)<=bedPtsZ(2)
    bedRangeZ = [bedPtsZ(1)-abs(voxDimZ)/2 bedPtsZ(2)+abs(voxDimZ)/2];
else
    bedRangeZ = [bedPtsZ(1)+abs(voxDimZ)/2 bedPtsZ(2)-abs(voxDimZ)/2];
end

end

function [IOP, IPP, imSz1, imSz2, voxDim1, voxDim2, varargout] = parseInputs(varargin)

if nargin>1
    IOP = varargin{1};
    IPP = varargin{2};
    imSz1 = varargin{3};
    imSz2 = varargin{4};
    voxDim1 = varargin{5};
    voxDim2 = varargin{6};
    if nargin>6
        voxDim3 = varargin{7};
    end
    if ~exist('IOP','var')
        error('Image Orientation Patient is not defined')
    end
    if ~exist('IPP','var')
        error('Image Position Patient is not defined')
    end
    if ~exist('imSz1','var')
        error('Image matrix size 1 is not defined')
    end
    if ~exist('imSz2','var')
        error('Image matrix size 2 is not defined')
    end
    if ~exist('voxDim1','var')
        error('Pixel size 1 is not defined')
    end
    if ~exist('voxDim2','var')
        error('Pixel size 2 is not defined')
    end
elseif isstruct(varargin{1})
    S = varargin{1};
    if isfield(S,'IOP')
        IOP = S.IOP;
    else
        error('Image Orientation Patient is not defined in input structure')
    end
    if isfield(S,'IPP')
        IPP = S.IPP;
    else
        error('Image Position Patient is not defined in input structure')
    end
    if isfield(S,'volumes')
        imSz1 = size(S.volumes,1);
        imSz2 = size(S.volumes,2);
    elseif isfield(S,'volume')
        imSz1 = size(S.volume,1);
        imSz2 = size(S.volume,2);
    else
        error('Image volume is not defined in input structure')
    end
    if isfield(S,'voxDim1')
        voxDim1 = S.voxDim1;
    elseif isfield(S,'bedRange1')
        voxDim1 = diff(S.bedRange1)/imSz1;
    else
        error('Cannot get voxel dimension 1 from input structure')
    end
    if isfield(S,'voxDim2')
        voxDim2 = S.voxDim2;
    elseif isfield(S,'bedRange2')
        voxDim2 = diff(S.bedRange2)/imSz2;
    else
        error('Cannot get voxel dimension 2 from input structure')
    end
    if isfield(S,'voxDim3')
        voxDim3 = S.voxDim3;
    elseif isfield(S,'bedRange3')
        try
            voxDim3 = diff(S.bedRange3)/size(S.volume,3);
        catch
            voxDim3 = diff(S.bedRange3)/size(S.volumes,3);
        end
    else
        error('Cannot get voxel dimension 3 from input structure')
    end
end

if exist('voxDim3','var'), varargout{1} = voxDim3; end

end
