function [locX, locY, locZ] = getVoxelLocs(varargin)

%this function calculates the complete set of voxel locations
%input is one corner column, like dicom field 'ImagePositionPatient'
%NOTE: format should be like that of 'readImages' structure,
%   NOT dicom, so loc1 counts X axis and loc2 counts Y axis

[IOP, IPP, imSz1, imSz2, imSz3, voxDim1, voxDim2] = parseImputs(varargin{:});

xr = IOP(1);
yr = IOP(2);
zr = IOP(3);
xc = IOP(4);
yc = IOP(5);
zc = IOP(6);
xs = yr*zc-yc*zr;
ys = zr*xc-zc*xr;
zs = xr*yc-xc*yr;
rot = [xr xc xs; yr yc ys; zr zc zs];

%construct volumes of absolute voxel positions from old IPP
[temp2,temp1] = meshgrid(0:imSz2-1,0:imSz1-1);
%these offsets will be the same for each slice
offsets = rot(:,1) * temp1(:)' * voxDim1 + rot(:,2) * temp2(:)' * voxDim2;
offsetX = reshape(offsets(1,:),imSz1,imSz2);
offsetY = reshape(offsets(2,:),imSz1,imSz2);
offsetZ = reshape(offsets(3,:),imSz1,imSz2);

%voxel locations are the base IPP, duplicated for each slice, plus offsets
try
    locX = reshape(repmat(IPP(1,:),imSz1*imSz2,1),imSz1,imSz2,imSz3)+...
        repmat(offsetX,1,1,imSz3);
    locY = reshape(repmat(IPP(2,:),imSz1*imSz2,1),imSz1,imSz2,imSz3)+...
        repmat(offsetY,1,1,imSz3);
    locZ = reshape(repmat(IPP(3,:),imSz1*imSz2,1),imSz1,imSz2,imSz3)+...
        repmat(offsetZ,1,1,imSz3);
catch
    tempX = zeros([size(offsetX) imSz3]);
    tempY = zeros([size(offsetY) imSz3]);
    tempZ = zeros([size(offsetZ) imSz3]);
    for n=1:imSz3
        tempX(:,:,n) = offsetX;
        tempY(:,:,n) = offsetY;
        tempZ(:,:,n) = offsetZ;
    end    
    locX = reshape(repmat(IPP(1,:),imSz1*imSz2,1),imSz1,imSz2,imSz3)+tempX;
    locY = reshape(repmat(IPP(2,:),imSz1*imSz2,1),imSz1,imSz2,imSz3)+tempY;
    locZ = reshape(repmat(IPP(3,:),imSz1*imSz2,1),imSz1,imSz2,imSz3)+tempZ;
end

end

function [IOP, IPP, imSz1, imSz2, imSz3, voxDim1, voxDim2] = parseImputs(varargin)

if nargin==7
    IOP = varargin{1};
    IPP = varargin{2};
    imSz1 = varargin{3};
    imSz2 = varargin{4};
    imSz3 = varargin{5};
    voxDim1 = varargin{6};
    voxDim2 = varargin{7};
elseif isstruct(varargin{1})
    IOP = varargin{1}.IOP;
    IPP = varargin{1}.IPP;
    if isfield(varargin{1},'vol') %precedence goes to 'vol'
        imSz1 = size(varargin{1}.vol,1);
        imSz2 = size(varargin{1}.vol,2);
        imSz3 = size(varargin{1}.vol,3);
    elseif isfield(varargin{1},'volume')
        imSz1 = size(varargin{1}.volume,1);
        imSz2 = size(varargin{1}.volume,2);
        imSz3 = size(varargin{1}.volume,3);
    elseif isfield(varargin{1},'volumes')
        imSz1 = size(varargin{1}.volumes,1);
        imSz2 = size(varargin{1}.volumes,2);
        imSz3 = size(varargin{1}.volumes,3);
    end
    if isfield(varargin{1},'voxDim1')
        voxDim1 = varargin{1}.voxDim1;
    else
        voxDim1 = diff(varargin{1}.bedRange1)/imSz1;
    end
    if isfield(varargin{1},'voxDim2')
        voxDim2 = varargin{1}.voxDim2;
    else
        voxDim2 = diff(varargin{1}.bedRange2)/imSz2;
    end
else
    error('Incorrect inputs.')
end

%assure single precision operation
voxDim1 = single(voxDim1);
voxDim2 = single(voxDim2);

end