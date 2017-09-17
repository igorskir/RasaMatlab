function smVol = GaussianSmooth(vol, FWHM)

%input is original volume and smoothing kernel FWHM (in voxel units)

%only 2D or 3D image data are supported
dim = 3;
if size(vol,3)==1, dim = 2; end

if nargin==1
    FWHM = 5; %default 5mm
end

if numel(FWHM)==1
    FWHM = repmat(FWHM,1,dim);
elseif numel(FWHM)>=dim
    FWHM = FWHM(1:dim);
else
    error(['Number of elements in smoothing size should be 1 or ' num2str(dim)])
end

FWHM = abs(FWHM(:)');

if ~all(FWHM==0)
    %gaussian kernel
    sig = FWHM/(2*sqrt(2*0.693147180559945));
    %radius width of kernal is 2X standard deviation (for each tail, ie. 4X total)
    siz = round(2*FWHM);
    if dim==2
        [x,y] = ndgrid(-siz(1):siz(1),-siz(2):siz(2));
        h = exp(-(x.*x/2/sig(1)^2 + y.*y/2/sig(2)^2));
    elseif dim==3
        [x,y,z] = ndgrid(-siz(1):siz(1),-siz(2):siz(2),-siz(3):siz(3));
        h = exp(-(x.*x/2/sig(1)^2 + y.*y/2/sig(2)^2 + z.*z/2/sig(3)^2));
    end
    h = h/sum(h(:));
    %     %convolve
    %     smVol = imfilter(vol, h, 'replicate');
    %copy methods from function 'smooth3'
    padSz = siz(:)';
    smVol = convn(padReplicate(vol,padSz),h,'valid');
else
    smVol = vol;
end

end

function out = padReplicate(data, padSz)

numDims = length(padSz);
idx = cell(numDims,1);
for n=1:numDims
    M = size(data,n);
    onesVector = ones(1,padSz(n));
    idx{n} = [onesVector 1:M M*onesVector];
end
out = data(idx{:});

end

