function [bedRange1, bedRange2, bedRange3, varargout] = getVolumeBedRanges(varargin)

%this function gets the bed coverage from IOP, IPP, image matrix size, and pixel spacing
%the output is in IMAGE coordinates, not SCANNER coordinates unless IOP==[1; 0; 0; 0; 1; 0]
%NOTE: voxDim1, i.e. PixelSpacing(1), is for rows (Y-axis when IOP==[1; 0; 0 0; 1; 0])
%      voxDim2, i.e. PixelSpacing(2), is for columns (X-axis when IOP==[1; 0; 0; 0; 1; 0])

[IOP, IPP, imSz1, imSz2, imSz3, voxDim1, voxDim2, voxDim3, locX, locY, locZ] = parseInputs(varargin{:});

%rotation matrix to put image coordinates into absolute scanner coordinates
vect1 = IOP(1:3);
vect2 = IOP(4:6);
vect3 = cross(vect1, vect2);
%construct volumes of absolute voxel positions
if isempty(locX) || isempty(locY) || isempty(locZ)
    [locX, locY, locZ] = getVoxelLocs(IOP,IPP,imSz1,imSz2,imSz3,voxDim1,voxDim2);
end
%calculate bed ranges (along image axes, which isn't scanner axes unless IOP==[1;0;0;0;1;0])
cntr1 = round(imSz1/2);
cntr2 = round(imSz2/2);
cntr3 = round(imSz3/2);
%bed starts and ends along new axes
bedPts1 = vect1' * [locX(1,cntr2,cntr3) locY(1,cntr2,cntr3) locZ(1,cntr2,cntr3); ...
    locX(end,cntr2,cntr3) locY(end,cntr2,cntr3) locZ(end,cntr2,cntr3)]';
bedPts2 = vect2' * [locX(cntr1,1,cntr3) locY(cntr1,1,cntr3) locZ(cntr1,1,cntr3); ...
    locX(cntr1,end,cntr3) locY(cntr1,end,cntr3) locZ(cntr1,end,cntr3)]';
bedPts3 = vect3' * [locX(cntr1,cntr2,1) locY(cntr1,cntr2,1) locZ(cntr1,cntr2,1); ...
    locX(cntr1,cntr2,end) locY(cntr1,cntr2,end) locZ(cntr1,cntr2,end)]';

%define bed ranges, use 'abs' method
if bedPts1(1)<=bedPts1(2)
    bedRange1 = [bedPts1(1)-abs(voxDim1)/2 bedPts1(2)+abs(voxDim1)/2];
else
    bedRange1 = [bedPts1(1)+abs(voxDim1)/2 bedPts1(2)-abs(voxDim1)/2];
end
if bedPts2(1)<=bedPts2(2)
    bedRange2 = [bedPts2(1)-abs(voxDim2)/2 bedPts2(2)+abs(voxDim2)/2];
else
    bedRange2 = [bedPts2(1)+abs(voxDim2)/2 bedPts2(2)-abs(voxDim2)/2];
end
if bedPts3(1)<=bedPts3(2)
    bedRange3 = [bedPts3(1)-abs(voxDim3)/2 bedPts3(2)+abs(voxDim3)/2];
else
    bedRange3 = [bedPts3(1)+abs(voxDim3)/2 bedPts3(2)-abs(voxDim3)/2];
end

%optional output full voxel maps
varargout{1} = locX;
varargout{2} = locY;
varargout{3} = locZ;

end

function [IOP, IPP, imSz1, imSz2, imSz3, voxDim1, voxDim2, voxDim3, locX, locY, locZ] = parseInputs(varargin)

%these will be changed if found in structure input
locX = [];
locY = [];
locZ = [];

if nargin>1
    IOP = varargin{1};
    IPP = varargin{2};
    imSz1 = varargin{3};
    imSz2 = varargin{4};
    imSz3 = size(IPP,2);
    voxDim1 = varargin{5};
    voxDim2 = varargin{6};
    if nargin>6
        voxDim3 = varargin{7};
    end
    if ~exist('IOP','var')
        error('Image Orientation Patient is not defined.')
    end
    if ~exist('IPP','var')
        error('Image Position Patient is not defined.')
    end
    if ~exist('imSz1','var')
        error('Image matrix size 1 is not defined.')
    end
    if ~exist('imSz2','var')
        error('Image matrix size 2 is not defined.')
    end
    if ~exist('voxDim1','var')
        error('Pixel size 1 is not defined.')
    end
    if ~exist('voxDim2','var')
        error('Pixel size 2 is not defined.')
    end
elseif isstruct(varargin{1})
    S = varargin{1};
    if isfield(S,'IOP')
        IOP = S.IOP;
    else
        error('Image Orientation Patient is not defined in input structure.')
    end
    if isfield(S,'IPP')
        IPP = S.IPP;
    else
        error('Image Position Patient is not defined in input structure.')
    end
    if isfield(S,'volumes')
        imSz1 = size(S.volumes,1);
        imSz2 = size(S.volumes,2);
    elseif isfield(S,'volume')
        imSz1 = size(S.volume,1);
        imSz2 = size(S.volume,2);
    elseif isfield(S,'vol')
        imSz1 = size(S.vol,1);
        imSz2 = size(S.vol,2);
    else
        error('Image volume is not defined in input structure.')
    end
    imSz3 = size(IPP,2);
    if isfield(S,'bedRange1')
        voxDim1 = diff(S.bedRange1)/imSz1;
    elseif isfield(S,'voxDim1')
        voxDim1 = S.voxDim1;
    else
        error('Cannot get voxel dimension 1 from input structure.')
    end
    if isfield(S,'bedRange2')
        voxDim2 = diff(S.bedRange2)/imSz2;
    elseif isfield(S,'voxDim2')
        voxDim2 = S.voxDim2;
    else
        error('Cannot get voxel dimension 2 from input structure.')
    end
    if isfield(S,'bedRange3')
        if isfield(S,'volume')
            voxDim3 = diff(S.bedRange3)/size(S.volume,3);
        elseif isfield(S,'volumes')
            voxDim3 = diff(S.bedRange3)/size(S.volumes,3);
        elseif isfield(S,'vol')
            voxDim3 = diff(S.bedRange3)/size(S.vol,3);
        end
    elseif isfield(S,'voxDim3')
        voxDim3 = S.voxDim3;
    else
        error('Cannot get voxel dimension 3 from input structure.')
    end
    if isfield(S,'locX')
        locX = S.locX;
    end
    if isfield(S,'locY')
        locY = S.locY;
    end
    if isfield(S,'locZ')
        locZ = S.locZ;
    end
end

if ~exist('voxDim3','var')
    if size(IPP,2)>1
        %most common slice interval, in case missing slices
        voxDim3 = mode(sqrt(sum(diff(IPP,[],2).^2)));
    else %single slice
        error('Must input thickness of a single slice.')
    end
end

end
