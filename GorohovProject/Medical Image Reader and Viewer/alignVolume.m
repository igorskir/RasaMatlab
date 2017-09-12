function out = alignVolume(out, varargin)

%this function aligns volume to (cosine) vector (i.e. [ImageOrientationPatient]')
%it may involve reslicing, flipping, or both
%input is structure from 'readImages' function
%varargin can be [IOP] or 'align'

[volume, IOP, voxDim1, voxDim2, voxDim3] = parseInput(out);

align = false; %align to scanner coordinate axes, including flipping?
verbose = false; %display update statements
if ~isempty(varargin)
    for n=1:length(varargin)
        if (isnumeric(varargin{n}) && ...
                (all(size(varargin{n})==[6 1]) || all(size(varargin{n})==[1 6])))
            newIOP = varargin{:}(:);
            %assure normalized vectors
            newIOP(1:3) = newIOP(1:3)/norm(newIOP(1:3));
            newIOP(4:6) = newIOP(4:6)/norm(newIOP(4:6));
            %check if valid IOP vector
            if any(isnan(newIOP)) || (round(newIOP(1:3)'*newIOP(4:6)*1e6)/1e6)~=0
                error('Invalid patient orientation input vector');
            end
        end
        if ischar(varargin{n})
            if strcmpi(varargin{n}(1),'a')
                align = true;
            elseif strcmpi(varargin{n}(1),'v')
                verbose = true;
            end
        end
    end
end

%if no orientation vector input, assume to align to scanner coordinate axes
if ~exist('newIOP','var')
    newIOP = [1;0;0;0;1;0];
    align = 1;
end

%reslice volume if not aligned to new IOP
if ~all(round(out.IOP*1e6)/1e6==round(newIOP*1e6)/1e6)
    if verbose, disp('   Reslicing...'); end
    temp = resliceVolume(out,newIOP);
    volume = temp.volume;
    voxDim1 = temp.voxDim1;
    voxDim2 = temp.voxDim2;
    voxDim3 = temp.voxDim3;
    IOP = temp.IOP;
    imSz1 = size(volume,1);
    imSz2 = size(volume,2);
    %bed ranges are along new IOP now
    bedRange1 = temp.bedRange1;
    bedRange2 = temp.bedRange2;
    bedRange3 = temp.bedRange3;
    out.locX = temp.locX;
    out.locY = temp.locY;
    out.locZ = temp.locZ;
else
    imSz1 = size(volume,1);
    imSz2 = size(volume,2);
    bedRange1 = out.bedRange1;
    bedRange2 = out.bedRange2;
    bedRange3 = out.bedRange3;
    if ~(isfield(out,'locX') && isfield(out,'locY') && isfield(out,'locZ'))
        [out.locX, out.locY, out.locZ] = getVoxelLocs(out);
    end
end
imSz3 = size(volume,3);

if align
    %now flip if needed, accounting for IOP
    cntr1 = round(imSz1/2);
    cntr2 = round(imSz2/2);
    cntr3 = round(imSz3/2);
    %rotation vectors
    vect1 = IOP(1:3);
    vect2 = IOP(4:6);
    vect3 = cross(vect1, vect2);
    if bedRange1(1)>bedRange1(2)
        if verbose, disp('   Flipping X dimension...'); end
        try
            volume = flip(volume,1);
        catch
            volume = flipdim(volume,1);
        end
        %get new bed ranges
        try
            out.locX = flip(out.locX,1);
            out.locY = flip(out.locY,1);
            out.locZ = flip(out.locZ,1);
        catch
            out.locX = flipdim(out.locX,1);
            out.locY = flipdim(out.locY,1);
            out.locZ = flipdim(out.locZ,1);
        end
        bedPts1 = vect1' * [out.locX(1,cntr2,cntr3) out.locY(1,cntr2,cntr3) out.locZ(1,cntr2,cntr3); ...
            out.locX(end,cntr2,cntr3) out.locY(end,cntr2,cntr3) out.locZ(end,cntr2,cntr3)]';
        if bedPts1(1)<=bedPts1(2)
            bedRange1 = [bedPts1(1)-abs(voxDim1)/2 bedPts1(2)+abs(voxDim1)/2];
        else
            bedRange1 = [bedPts1(1)+abs(voxDim1)/2 bedPts1(2)-abs(voxDim1)/2];
        end
        voxDim1 = -voxDim1;
    end
    if bedRange2(1)>bedRange2(2)
        if verbose, disp('   Flipping Y dimension...'); end
        try
            volume = flip(volume,2);
        catch
            volume = flipdim(volume,2);
        end
        %get new bed ranges
        try
            out.locX = flip(out.locX,2);
            out.locY = flip(out.locY,2);
            out.locZ = flip(out.locZ,2);
        catch
            out.locX = flipdim(out.locX,2);
            out.locY = flipdim(out.locY,2);
            out.locZ = flipdim(out.locZ,2);
        end
        bedPts2 = vect2' * [out.locX(cntr1,1,cntr3) out.locY(cntr1,1,cntr3) out.locZ(cntr1,1,cntr3); ...
            out.locX(cntr1,end,cntr3) out.locY(cntr1,end,cntr3) out.locZ(cntr1,end,cntr3)]';
        if bedPts2(1)<=bedPts2(2)
            bedRange2 = [bedPts2(1)-abs(voxDim2)/2 bedPts2(2)+abs(voxDim2)/2];
        else
            bedRange2 = [bedPts2(1)+abs(voxDim2)/2 bedPts2(2)-abs(voxDim2)/2];
        end
        voxDim2 = -voxDim2;
    end
    if bedRange3(1)>bedRange3(2)
        if verbose, disp('   Flipping Z dimension...'); end
        try
            volume = flip(volume,3);
        catch
            volume = flipdim(volume,3);
        end
        %get new bed ranges
        try
            out.locX = flip(out.locX,3);
            out.locY = flip(out.locY,3);
            out.locZ = flip(out.locZ,3);
        catch
            out.locX = flipdim(out.locX,3);
            out.locY = flipdim(out.locY,3);
            out.locZ = flipdim(out.locZ,3);
        end
        bedPts3 = vect3' * [out.locX(cntr1,cntr2,1) out.locY(cntr1,cntr2,1) out.locZ(cntr1,cntr2,1); ...
            out.locX(cntr1,cntr2,end) out.locY(cntr1,cntr2,end) out.locZ(cntr1,cntr2,end)]';
        if bedPts3(1)<=bedPts3(2)
            bedRange3 = [bedPts3(1)-abs(voxDim3)/2 bedPts3(2)+abs(voxDim3)/2];
        else
            bedRange3 = [bedPts3(1)+abs(voxDim3)/2 bedPts3(2)-abs(voxDim3)/2];
        end
        voxDim3 = -voxDim3;
    end
end

if isfield(out,'volume')
    out.volume = volume;
elseif isfield(out,'volumes')
    out.volumes = volume;
end
out.imSz1 = size(volume,1);
out.imSz2 = size(volume,2);
out.imSz3 = size(volume,3);
out.voxDim1 = voxDim1;
out.voxDim2 = voxDim2;
out.voxDim3 = voxDim3;
out.bedRange1 = bedRange1;
out.bedRange2 = bedRange2;
out.bedRange3 = bedRange3;
out.IOP = IOP;
%redefine IPP
out.IPP = [squeeze(out.locX(1,1,1:end))'; squeeze(out.locY(1,1,1:end))'; squeeze(out.locZ(1,1,1:end))'];

%save memory
out.locX = single(out.locX);
out.locY = single(out.locY);
out.locZ = single(out.locZ);

if isfield(out,'patientOrientation')
    out = rmfield(out,'patientOrientation');
end

end

function [volume, IOP, voxDim1, voxDim2, voxDim3] = parseInput(in)

if isfield(in,'volumes') || isfield(in,'volume')
    if isfield(in,'volumes')
        volume = in.volumes;
    else
        volume = in.volume;
    end
else
    error('Image volume must be defined')
end
if isfield(in,'IOP')
    IOP = in.IOP;
else
    error('Image orientation relative to scanner must be defined')
end
if isfield(in,'voxDim1')
    voxDim1 = in.voxDim1;
elseif isfield(in,'bedRange1')
    voxDim1 = diff(in.bedRange1)/size(volume,1);
else
    error('Cannot get voxel dimension 1 from input structure')
end
if isfield(in,'voxDim2')
    voxDim2 = in.voxDim2;
elseif isfield(in,'bedRange2')
    voxDim2 = diff(in.bedRange2)/size(volume,2);
else
    error('Cannot get voxel dimension 2 from input structure')
end
if isfield(in,'voxDim3')
    voxDim3 = in.voxDim3;
elseif isfield(in,'bedRange3')
    voxDim3 = diff(in.bedRange3)/size(volume,3);
else
    error('Cannot get voxel dimension 3 from input structure')
end

if size(volume,3)==1
    warning('CAREFUL: Strange things happen when rotating slices')
end

end
