function out = coregisterLocs(in1, in2, varargin)

%this function returns voxel location matrices for the coregistration
%the volumes will be aligned to the orientation of the first input volume
%only coregistered output points are calculated here,
%   input points come from 'rotateInterpolationPoints' function

trim = false; %show all volume space or trim to show only shared space?
align = false; %align to scanner coordinate axes?
full = false; %resample image matrices, regardless of how large?

if ~isempty(varargin)
    for n=1:length(varargin)
        if ischar(varargin{n}) && strcmpi(varargin{n}(1),'t')
            trim = true;
        end
        if ischar(varargin{n}) && strcmpi(varargin{n}(1),'a')
            align = true;
        end
        if ischar(varargin{n}) && strcmpi(varargin{n}(1),'f')
            full = true;
        end
    end
end

%throw error if one set is 2D and the other is 3D
if (size(in1.volumes,3)==1 && size(in2.volumes,3)>1) || ...
        (size(in1.volumes,3)>1 && size(in2.volumes,3)==1)
    disp('Inputs have different dimensions and should not be coregistered')
    return
end
%2D data need additional checks
if (size(in1.volumes,3)==1 && size(in2.volumes,3)==1)
    if~all(in1.IOP==in2.IOP)
        %throw error if both sets are 2D and have different orientations
        disp('Input slices have different orientations and cannot be coregistered')
        return
    end
    if align
        warning('CAREFUL: Strange things happen when rotating slices')
    end
end

%volumes needed for input to intperpolation function
if isfield(in1,'volume')
    out.vol_1 = in1.volume;
elseif isfield(in1,'volumes')
    out.vol_1 = in1.volumes;
end
if isfield(in2,'volume')
    out.vol_2 = in2.volume;
elseif isfield(in2,'volumes')
    out.vol_2 = in2.volumes;
end

oldImStempZ1(1) = size(out.vol_1,1);
oldImStempZ1(end) = size(out.vol_1,2);
oldImSz3_1 = size(out.vol_1,3);
oldImStempZ2(1) = size(out.vol_2,1);
oldImStempZ2(end) = size(out.vol_2,2);
oldImSz3_2 = size(out.vol_2,3);

if ~(isfield(in1,'locX') && isfield(in1,'locY') && isfield(in1,'locZ'))
    [in1.locX,in1.locY,in1.locZ] = getVoxelLocs(in1);
end

%get old rotations of each volume
oldRot_1 = [in1.IOP(1:3) in1.IOP(4:6) cross(in1.IOP(1:3),in1.IOP(4:6))];
oldRot_2 = [in2.IOP(1:3) in2.IOP(4:6) cross(in2.IOP(1:3),in2.IOP(4:6))];

if align
    out.IOP = [1;0;0;0;1;0];
    newRot = eye(3);
else
    %align both to frame of first volume
    out.IOP = in1.IOP;
    newRot = oldRot_1;
end

if all(out.IOP==in1.IOP) && all(in1.IOP==in2.IOP) && (oldImStempZ1(1)==oldImStempZ2(1) && oldImStempZ1(end)==oldImStempZ2(end) && oldImSz3_1==oldImSz3_2) && ...
        (all(in1.bedRange1==in2.bedRange1) && all(in1.bedRange2==in2.bedRange2) && all(in1.bedRange3==in2.bedRange3))
    %volumes are already have exactly the same sampling
    [out.oldInterpX_1,out.oldInterpY_1,out.oldInterpZ_1] = getRotatedInterpolationGrid(in1.locX,in1.locY,in1.locZ,in1.IOP,oldRot_1');
    %assure NDGRID format
    [out.oldInterpY_1,out.oldInterpX_1,out.oldInterpZ_1] = ...
        meshgrid(squeeze(out.oldInterpY_1(1,1:end,1)),squeeze(out.oldInterpX_1(1:end,1,1)),squeeze(out.oldInterpZ_1(1,1,1:end)));
    out.oldInterpX_2 = out.oldInterpX_1;
    out.oldInterpY_2 = out.oldInterpY_1;
    out.oldInterpZ_2 = out.oldInterpZ_1;
    out.newInterpX_1 = out.oldInterpX_1;
    out.newInterpY_1 = out.oldInterpY_1;
    out.newInterpZ_1 = out.oldInterpZ_1;
    out.newInterpX_2 = out.oldInterpX_1;
    out.newInterpY_2 = out.oldInterpY_1;
    out.newInterpZ_2 = out.oldInterpZ_1;
    out.locX = in1.locX;
    out.locY = in1.locY;
    out.locZ = in1.locZ;
    out.IOP = in1.IOP;
    out.IPP = in1.IPP;
    out.voxDim1 = diff(in1.bedRange1)/size(out.locX,1);
    out.voxDim2 = diff(in1.bedRange2)/size(out.locX,2);
    out.voxDim3 = diff(in1.bedRange3)/size(out.locX,3);
    
    %assure correct relative orientation if 'align' flag set
    if align
        if out.voxDim1<0
            try
                out.locX = flip(out.locX,1);
                out.locY = flip(out.locY,1);
                out.locZ = flip(out.locZ,1);
                out.newInterpX_1 = flip(out.newInterpX_1,1);
                out.newInterpY_1 = flip(out.newInterpY_1,1);
                out.newInterpZ_1 = flip(out.newInterpZ_1,1);
                out.newInterpX_2 = flip(out.newInterpX_2,1);
                out.newInterpY_2 = flip(out.newInterpY_2,1);
                out.newInterpZ_2 = flip(out.newInterpZ_2,1);
            catch
                out.locX = flipdim(out.locX,1);
                out.locY = flipdim(out.locY,1);
                out.locZ = flipdim(out.locZ,1);
                out.newInterpX_1 = flipdim(out.newInterpX_1,1);
                out.newInterpY_1 = flipdim(out.newInterpY_1,1);
                out.newInterpZ_1 = flipdim(out.newInterpZ_1,1);
                out.newInterpX_2 = flipdim(out.newInterpX_2,1);
                out.newInterpY_2 = flipdim(out.newInterpY_2,1);
                out.newInterpZ_2 = flipdim(out.newInterpZ_2,1);
            end
            out.voxDim1 = -out.voxDim1;
        end
        if out.voxDim2<0
            try
                out.locX = flip(out.locX,2);
                out.locY = flip(out.locY,2);
                out.locZ = flip(out.locZ,2);
                out.newInterpX_1 = flip(out.newInterpX_1,2);
                out.newInterpY_1 = flip(out.newInterpY_1,2);
                out.newInterpZ_1 = flip(out.newInterpZ_1,2);
                out.newInterpX_2 = flip(out.newInterpX_2,2);
                out.newInterpY_2 = flip(out.newInterpY_2,2);
                out.newInterpZ_2 = flip(out.newInterpZ_2,2);
            catch
                out.locX = flipdim(out.locX,2);
                out.locY = flipdim(out.locY,2);
                out.locZ = flipdim(out.locZ,2);
                out.newInterpX_1 = flipdim(out.newInterpX_1,2);
                out.newInterpY_1 = flipdim(out.newInterpY_1,2);
                out.newInterpZ_1 = flipdim(out.newInterpZ_1,2);
                out.newInterpX_2 = flipdim(out.newInterpX_2,2);
                out.newInterpY_2 = flipdim(out.newInterpY_2,2);
                out.newInterpZ_2 = flipdim(out.newInterpZ_2,2);
            end
            out.voxDim2 = -out.voxDim2;
        end
        if out.voxDim3<0
            try
                out.locX = flip(out.locX,3);
                out.locY = flip(out.locY,3);
                out.locZ = flip(out.locZ,3);
                out.newInterpX_1 = flip(out.newInterpX_1,3);
                out.newInterpY_1 = flip(out.newInterpY_1,3);
                out.newInterpZ_1 = flip(out.newInterpZ_1,3);
                out.newInterpX_2 = flip(out.newInterpX_2,3);
                out.newInterpY_2 = flip(out.newInterpY_2,3);
                out.newInterpZ_2 = flip(out.newInterpZ_2,3);
            catch
                out.locX = flipdim(out.locX,3);
                out.locY = flipdim(out.locY,3);
                out.locZ = flipdim(out.yocZ,3);
                out.newInterpX_1 = flipdim(out.newInterpX_1,3);
                out.newInterpY_1 = flipdim(out.newInterpY_1,3);
                out.newInterpZ_1 = flipdim(out.newInterpZ_1,3);
                out.newInterpX_2 = flipdim(out.newInterpX_2,3);
                out.newInterpY_2 = flipdim(out.newInterpY_2,3);
                out.newInterpZ_2 = flipdim(out.newInterpZ_2,3);
            end
            out.voxDim3 = -out.voxDim3;
        end
    end
    out.IPP = [squeeze(out.locX(1,1,1:end))'; ...
        squeeze(out.locY(1,1,1:end))'; ...
        squeeze(out.locZ(1,1,1:end))'];
    return
end

%won't use new individual sampling points so pass 'coreg' flag
temp1 = rotateInterpolationPoints(in1,newRot,'coreg');
temp2 = rotateInterpolationPoints(in2,newRot,'coreg');
%'old' and 'new' output locations are not correct but are in the relative same space
clear in1 in2

%assure correct relative orientation if 'align' flag set
if align
    if temp1.newVoxDim1<0
        try
            temp1.newLocX = flip(temp1.newLocX,1);
            temp1.newLocY = flip(temp1.newLocY,1);
            temp1.newLocZ = flip(temp1.newLocZ,1);
        catch
            temp1.newLocX = flipdim(temp1.newLocX,1);
            temp1.newLocY = flipdim(temp1.newLocY,1);
            temp1.newLocZ = flipdim(temp1.newLocZ,1);
        end
        temp1.newVoxDim1 = -temp1.newVoxDim1;
    end
    if temp1.newVoxDim2<0
        try
            temp1.newLocX = flip(temp1.newLocX,2);
            temp1.newLocY = flip(temp1.newLocY,2);
            temp1.newLocZ = flip(temp1.newLocZ,2);
        catch
            temp1.newLocX = flipdim(temp1.newLocX,2);
            temp1.newLocY = flipdim(temp1.newLocY,2);
            temp1.newLocZ = flipdim(temp1.newLocZ,2);
        end
        temp1.newVoxDim2 = -temp1.newVoxDim2;
    end
    if temp1.newVoxDim3<0
        try
            temp1.newLocX = flip(temp1.newLocX,3);
            temp1.newLocY = flip(temp1.newLocY,3);
            temp1.newLocZ = flip(temp1.newLocZ,3);
        catch
            temp1.newLocX = flipdim(temp1.newLocX,3);
            temp1.newLocY = flipdim(temp1.newLocY,3);
            temp1.newLocZ = flipdim(temp1.newLocZ,3);
        end
        temp1.newVoxDim3 = -temp1.newVoxDim3;
    end
    if temp2.newVoxDim1<0
        try
            temp2.newLocX = flip(temp2.newLocX,1);
            temp2.newLocY = flip(temp2.newLocY,1);
            temp2.newLocZ = flip(temp2.newLocZ,1);
        catch
            temp2.newLocX = flipdim(temp2.newLocX,1);
            temp2.newLocY = flipdim(temp2.newLocY,1);
            temp2.newLocZ = flipdim(temp2.newLocZ,1);
        end
        temp2.newVoxDim1 = -temp2.newVoxDim1;
    end
    if temp2.newVoxDim2<0
        try
            temp2.newLocX = flip(temp2.newLocX,2);
            temp2.newLocY = flip(temp2.newLocY,2);
            temp2.newLocZ = flip(temp2.newLocZ,2);
        catch
            temp2.newLocX = flipdim(temp2.newLocX,2);
            temp2.newLocY = flipdim(temp2.newLocY,2);
            temp2.newLocZ = flipdim(temp2.newLocZ,2);
        end
        temp2.newVoxDim2 = -temp2.newVoxDim2;
    end
    if temp2.newVoxDim3<0
        try
            temp2.newLocX = flip(temp2.newLocX,3);
            temp2.newLocY = flip(temp2.newLocY,3);
            temp2.newLocZ = flip(temp2.newLocZ,3);
        catch
            temp2.newLocX = flipdim(temp2.newLocX,3);
            temp2.newLocY = flipdim(temp2.newLocY,3);
            temp2.newLocZ = flipdim(temp2.newLocZ,3);
        end
        temp2.newVoxDim3 = -temp2.newVoxDim3;
    end
else
    %assure relative orientation of 2nd volume matches that of 1st
    if temp1.newVoxDim1<0 && temp2.newVoxDim1>0 || ...
            temp1.newVoxDim1>0 && temp2.newVoxDim1<0
        try
            temp2.newLocX = flip(temp2.newLocX,1);
            temp2.newLocY = flip(temp2.newLocY,1);
            temp2.newLocZ = flip(temp2.newLocZ,1);
        catch
            temp2.newLocX = flipdim(temp2.newLocX,1);
            temp2.newLocY = flipdim(temp2.newLocY,1);
            temp2.newLocZ = flipdim(temp2.newLocZ,1);
        end
        temp2.newVoxDim1 = -temp2.newVoxDim1;
    end
    if temp1.newVoxDim2<0 && temp2.newVoxDim2>0 || ...
            temp1.newVoxDim2>0 && temp2.newVoxDim2<0
        try
            temp2.newLocX = flip(temp2.newLocX,2);
            temp2.newLocY = flip(temp2.newLocY,2);
            temp2.newLocZ = flip(temp2.newLocZ,2);
        catch
            temp2.newLocX = flipdim(temp2.newLocX,2);
            temp2.newLocY = flipdim(temp2.newLocY,2);
            temp2.newLocZ = flipdim(temp2.newLocZ,2);
        end
        temp2.newVoxDim2 = -temp2.newVoxDim2;
    end
    if (temp1.newVoxDim3<0 && temp2.newVoxDim3>0) || ...
            (temp1.newVoxDim3>0 && temp2.newVoxDim3<0)
        try
            temp2.newLocX = flip(temp2.newLocX,3);
            temp2.newLocY = flip(temp2.newLocY,3);
            temp2.newLocZ = flip(temp2.newLocZ,3);
        catch
            temp2.newLocX = flipdim(temp2.newLocX,3);
            temp2.newLocY = flipdim(temp2.newLocY,3);
            temp2.newLocZ = flipdim(temp2.newLocZ,3);
        end
        temp2.newVoxDim3 = -temp2.newVoxDim3;
    end
end

%check if there is overlap in the volumes, just need "aligned" corners
temp = newRot'*[temp1.newLocX(1,1,1) temp1.newLocY(1,1,1) temp1.newLocZ(1,1,1);...
    temp1.newLocX(end,1,1) temp1.newLocY(end,1,1) temp1.newLocZ(end,1,1);...
    temp1.newLocX(end,end,1) temp1.newLocY(end,end,1) temp1.newLocZ(end,end,1);...
    temp1.newLocX(1,end,1) temp1.newLocY(1,end,1) temp1.newLocZ(1,end,1);...
    temp1.newLocX(1,end,end) temp1.newLocY(1,end,end) temp1.newLocZ(1,end,end);...
    temp1.newLocX(1,1,end) temp1.newLocY(1,1,end) temp1.newLocZ(1,1,end);...
    temp1.newLocX(end,1,end) temp1.newLocY(end,1,end) temp1.newLocZ(end,1,end);...
    temp1.newLocX(end,end,end) temp1.newLocY(end,end,end) temp1.newLocZ(end,end,end)]';
tempX1 = temp(1,:);
tempY1 = temp(2,:);
tempZ1 = temp(3,:);
temp = newRot'*[temp2.newLocX(1,1,1) temp2.newLocY(1,1,1) temp2.newLocZ(1,1,1);...
    temp2.newLocX(end,1,1) temp2.newLocY(end,1,1) temp2.newLocZ(end,1,1);...
    temp2.newLocX(end,end,1) temp2.newLocY(end,end,1) temp2.newLocZ(end,end,1);...
    temp2.newLocX(1,end,1) temp2.newLocY(1,end,1) temp2.newLocZ(1,end,1);...
    temp2.newLocX(1,end,end) temp2.newLocY(1,end,end) temp2.newLocZ(1,end,end);...
    temp2.newLocX(1,1,end) temp2.newLocY(1,1,end) temp2.newLocZ(1,1,end);...
    temp2.newLocX(end,1,end) temp2.newLocY(end,1,end) temp2.newLocZ(end,1,end);...
    temp2.newLocX(end,end,end) temp2.newLocY(end,end,end) temp2.newLocZ(end,end,end)]';
tempX2 = temp(1,:);
tempY2 = temp(2,:);
tempZ2 = temp(3,:);

if min(tempX1(:))>max(tempX2(:)) || max(tempX1(:))<min(tempX2(:))
    if trim
        error('Volumes do not share the same space in the X dimension')
    else
        warning('Volumes do not share the same space in the X dimension')
    end
end
if min(tempY1(:))>max(tempY2(:)) || max(tempY1(:))<min(tempY2(:))
    if trim
        error('Volumes do not share the same space in the Y dimension')
    else
        warning('Volumes do not share the same space in the Y dimension')
    end
end
if min(tempZ1(:))>max(tempZ2(:)) || max(tempZ1(:))<min(tempZ2(:))
    if trim
        error('Volumes do not share the same space in the Z dimension')
    else
        warning('Volumes do not share the same space in the Z dimension')
    end
end
%check relative alignments (should never error)
if (tempX1(1)<tempX1(end) && tempX2(1)>tempX2(end)) || ...
        (tempX1(1)>tempX1(end) && tempX2(1)<tempX2(end)) || ...
        (tempY1(1)<tempY1(end) && tempY2(1)>tempY2(end)) || ...
        (tempY1(1)>tempY1(end) && tempY2(1)<tempY2(end)) || ...
        (tempZ1(1)<tempZ1(end) && tempZ2(1)>tempZ2(end)) || ...
        (tempZ1(1)>tempZ1(end) && tempZ2(1)<tempZ2(end))
    error('Input grids are not relatively aligned')
end

%old interpolation points come directly from 'rotateInterpolationPoints' function
out.oldInterpX_1 = temp1.oldInterpX;
out.oldInterpY_1 = temp1.oldInterpY;
out.oldInterpZ_1 = temp1.oldInterpZ;
out.oldInterpX_2 = temp2.oldInterpX;
out.oldInterpY_2 = temp2.oldInterpY;
out.oldInterpZ_2 = temp2.oldInterpZ;

%check for special case: same output voxel space, e.g. same volumes but different orientations
if all(size(temp1.newLocX)==size(temp2.newLocX))
    if temp1.newLocX(1,1,1)==temp2.newLocX(1,1,1) && ...
            temp1.newLocX(end,1,1)==temp2.newLocX(end,1,1) && ...
            temp1.newLocX(1,end,1)==temp2.newLocX(1,end,1) && ...
            temp1.newLocX(1,1,end)==temp2.newLocX(1,1,end) && ...
            temp1.newLocY(1,1,1)==temp2.newLocY(1,1,1) && ...
            temp1.newLocY(end,1,1)==temp2.newLocY(end,1,1) && ...
            temp1.newLocY(1,end,1)==temp2.newLocY(1,end,1) && ...
            temp1.newLocY(1,1,end)==temp2.newLocY(1,1,end) && ...
            temp1.newLocZ(1,1,1)==temp2.newLocZ(1,1,1) && ...
            temp1.newLocZ(end,1,1)==temp2.newLocZ(end,1,1) && ...
            temp1.newLocZ(1,end,1)==temp2.newLocZ(1,end,1) && ...
            temp1.newLocZ(1,1,end)==temp2.newLocZ(1,1,end)
        out.locX = temp1.newLocX;
        out.locY = temp1.newLocY;
        out.locZ = temp1.newLocZ;
        %calculate the new points differently, depending on if aligning to axes or not
        prec = 1e9;
        vectTest_1 = round((oldRot_1'*newRot)*prec)/prec;
        vectTest_2 = round((oldRot_2'*newRot)*prec)/prec;
        if all(vectTest_1(:)==0 | abs(vectTest_1(:))==1) && ...
                all(vectTest_2(:)==0 | abs(vectTest_2(:))==1)
            %this is just a flipping rotation and uses the original points
            [out.newInterpX_1, out.newInterpY_1, out.newInterpZ_1] = ...
                locFlip(temp1.oldInterpX,temp1.oldInterpY,temp1.oldInterpZ,vectTest_1,...
                [temp1.oldVoxDim1 temp1.oldVoxDim2 temp1.oldVoxDim3],...
                [temp1.newVoxDim1 temp1.newVoxDim2 temp1.newVoxDim3]);
            [out.newInterpX_2, out.newInterpY_2, out.newInterpZ_2] = ...
                locFlip(temp2.oldInterpX,temp2.oldInterpY,temp2.oldInterpZ,vectTest_2,...
                [temp2.oldVoxDim1 temp2.oldVoxDim2 temp2.oldVoxDim3],...
                [temp2.newVoxDim1 temp2.newVoxDim2 temp2.newVoxDim3]);
        else
            %this rotation requires newly calculated points
            [out.newInterpX_1, out.newInterpY_1, out.newInterpZ_1] = ...
                getRotatedInterpolationGrid(out.locX, out.locY, out.locZ, newRot(1:6), oldRot_1');
            [out.newInterpX_2, out.newInterpY_2, out.newInterpZ_2] = ...
                getRotatedInterpolationGrid(out.locX, out.locY, out.locZ, newRot(1:6), oldRot_2');
        end
        
        out.voxDim1 = temp1.newVoxDim1;
        out.voxDim2 = temp1.newVoxDim2;
        out.voxDim3 = temp1.newVoxDim3;
        
        out.IPP = [squeeze(out.locX(1,1,1:end))'; ...
            squeeze(out.locY(1,1,1:end))'; ...
            squeeze(out.locZ(1,1,1:end))'];
        return
    end
end

voxDim1_1 = temp1.newVoxDim1;
voxDim2_1 = temp1.newVoxDim2;
voxDim3_1 = temp1.newVoxDim3;
voxDim1_2 = temp2.newVoxDim1;
voxDim2_2 = temp2.newVoxDim2;
voxDim3_2 = temp2.newVoxDim3;

clear temp1 temp2

%find ranges of new combined interpolation (in input grid format)
if trim
    if tempX1(1)<tempX1(end)
        x1 = max([tempX1(1) tempX2(1)]);
        x2 = min([tempX1(end) tempX2(end)]);
    else
        x1 = min([tempX1(1) tempX2(1)]);
        x2 = max([tempX1(end) tempX2(end)]);
    end
    if tempY1(1)<tempY1(end)
        y1 = max([tempY1(1) tempY2(1)]);
        y2 = min([tempY1(end) tempY2(end)]);
    else
        y1 = min([tempY1(1) tempY2(1)]);
        y2 = max([tempY1(end) tempY2(end)]);
    end
    if tempZ1(1)<tempZ1(end)
        z1 = max([tempZ1(1) tempZ2(1)]);
        z2 = min([tempZ1(end) tempZ2(end)]);
    else
        z1 = min([tempZ1(1) tempZ2(1)]);
        z2 = max([tempZ1(end) tempZ2(end)]);
    end
else
    if tempX1(1)<tempX1(end)
        x1 = min([tempX1(1) tempX2(1)]);
        x2 = max([tempX1(end) tempX2(end)]);
    else
        x1 = max([tempX1(1) tempX2(1)]);
        x2 = min([tempX1(end) tempX2(end)]);
    end
    if tempY1(1)<tempY1(end)
        y1 = min([tempY1(1) tempY2(1)]);
        y2 = max([tempY1(end) tempY2(end)]);
    else
        y1 = max([tempY1(1) tempY2(1)]);
        y2 = min([tempY1(end) tempY2(end)]);
    end
    if tempZ1(1)<tempZ1(end)
        z1 = min([tempZ1(1) tempZ2(1)]);
        z2 = max([tempZ1(end) tempZ2(end)]);
    else
        z1 = max([tempZ1(1) tempZ2(1)]);
        z2 = min([tempZ1(end) tempZ2(end)]);
    end
end

%get voxel dimensions of combined interpolation
downSamp = 0; %default is to use highest pixel density
if downSamp
    %use coarsest voxels
    voxDim1abs = max(abs([voxDim1_1 voxDim1_2]));
    voxDim2abs = max(abs([voxDim2_1 voxDim2_2]));
    voxDim3abs = max(abs([voxDim3_1 voxDim3_2]));
else
    %use highest resolution
    voxDim1abs = min(abs([voxDim1_1 voxDim1_2]));
    voxDim2abs = min(abs([voxDim2_1 voxDim2_2]));
    voxDim3abs = min(abs([voxDim3_1 voxDim3_2]));
end

%generate output grids for rotation, use grid points from volume with chosen voxel spacing
%compare 2nd volume first to force potentially larger matrix (in the case of same volume)
%X
if voxDim1abs==abs(voxDim1_2)
    voxDim1 = voxDim1_2;
    newXpts = (tempX2(1):voxDim1:tempX2(end));
else
    voxDim1 = voxDim1_1;
    newXpts = (tempX1(1):voxDim1:tempX1(end));
end
%find points that lie within new bed range, include shared and padded ranges
if x1<x2
    temp = newXpts(newXpts>=x1&newXpts<=x2);
else
    temp = newXpts(newXpts<=x1&newXpts>=x2);
end
try
    pad1 = flip(temp(1)-voxDim1:-voxDim1:x1);
catch
    pad1 = fliplr(temp(1)-voxDim1:-voxDim1:x1);
end
pad2 = temp(end)+voxDim1:voxDim1:x2;
newXpts = [pad1 temp pad2];
%Y
if voxDim2abs==abs(voxDim2_2)
    voxDim2 = voxDim2_2;
    newYpts = (tempY2(1):voxDim2:tempY2(end));
else
    voxDim2 = voxDim2_1;
    newYpts = (tempY1(1):voxDim2:tempY1(end));
end
%find points that lie within new bed range, include shared and padded ranges
if y1<y2
    temp = newYpts(newYpts>=y1&newYpts<=y2);
else
    temp = newYpts(newYpts<=y1&newYpts>=y2);
end
try
    pad1 = flip(temp(1)-voxDim2:-voxDim2:y1);
catch
    pad1 = fliplr(temp(1)-voxDim2:-voxDim2:y1);
end
pad2 = temp(end)+voxDim2:voxDim2:y2;
newYpts = [pad1 temp pad2];
%Z
if voxDim3abs==abs(voxDim3_2)
    voxDim3 = voxDim3_2;
    newZpts = (tempZ2(1):voxDim3:tempZ2(end));
else
    voxDim3 = voxDim3_1;
    newZpts = (tempZ1(1):voxDim3:tempZ1(end));
end
%find points that lie within new bed range, include shared and padded ranges
if z1<z2
    temp = newZpts(newZpts>=z1&newZpts<=z2);
else
    temp = newZpts(newZpts<=z1&newZpts>=z2);
end
try
    pad1 = flip(temp(1)-voxDim3:-voxDim3:z1);
catch
    pad1 = fliplr(temp(1)-voxDim3:-voxDim3:z1);
end
pad2 = temp(end)+voxDim3:voxDim3:z2;
newZpts = [pad1 temp pad2];
clear temp

%avoid memory problems from working with huge matrices
if ~full
    flag = false;
    maxMatSz = 600;
    while numel(newXpts)>maxMatSz
        int = mean(unique(diff(newXpts)));
        temp = newXpts(1):int*2:newXpts(end);
        newXpts = temp;
        voxDim1 = voxDim1*2;
        voxDim1abs = voxDim1abs*2;
        flag = true;
    end
    while numel(newYpts)>maxMatSz
        int = mean(unique(diff(newYpts)));
        temp = newYpts(1):int*2:newYpts(end);
        newYpts = temp;
        voxDim2 = voxDim2*2;
        voxDim2abs = voxDim2abs*2;
        flag = true;
    end
    while numel(newZpts)>maxMatSz
        int = mean(unique(diff(newZpts)));
        temp = newZpts(1):int*2:newZpts(end);
        newZpts = temp;
        voxDim3 = voxDim3*2;
        voxDim3abs = voxDim3abs*2;
        flag = true;
    end
    if flag, disp('         using reduced output matrix sizes'); end
    clear temp int
end

%new volume sizes
newImSz1 = length(newXpts);
newImSz2 = length(newYpts);
newImSz3 = length(newZpts);

%offset to reshift points after rotation
%round to deal with precision errors when calculating isocenter from different orientations
isocntr = round(newRot*...
    [mean([newXpts(1) newXpts(end)]) ...
    mean([newYpts(1) newYpts(end)]) ...
    mean([newZpts(1) newZpts(end)])]'*1e2)/1e2; %precision to 0.01 mm
clear newXpts newYpts newZpts

%now create new grids of zero-centered sampling points, which are then reshifted and rotated for each volume
if x1<x2
    tempX = single(-voxDim1abs*(newImSz1-1)/2:voxDim1abs:voxDim1abs*(newImSz1-1)/2);
else
    tempX = single(voxDim1abs*(newImSz1-1)/2:-voxDim1abs:-voxDim1abs*(newImSz1-1)/2);
end
if y1<y2
    tempY = single(-voxDim2abs*(newImSz2-1)/2:voxDim2abs:voxDim2abs*(newImSz2-1)/2);
else
    tempY = single(voxDim2abs*(newImSz2-1)/2:-voxDim2abs:-voxDim2abs*(newImSz2-1)/2);
end
if z1<z2
    tempZ = single(-voxDim3abs*(newImSz3-1)/2:voxDim3abs:voxDim3abs*(newImSz3-1)/2);
else
    tempZ = single(voxDim3abs*(newImSz3-1)/2:-voxDim3abs:-voxDim3abs*(newImSz3-1)/2);
end


[newLocY,newLocX,newLocZ] = meshgrid(tempY,tempX,tempZ);
clear tempX tempY tempZ

%these are the newly rotated REAL locations in the output frame
[out.locX, out.locY, out.locZ] = ...
    getRotatedInterpolationGrid(newLocX, newLocY, newLocZ, [1;0;0;0;1;0], newRot);
%recenter
clear newLocX newLocY newLocZ
out.locX = out.locX + isocntr(1);
out.locY = out.locY + isocntr(2);
out.locZ = out.locZ + isocntr(3);

%get interpolation points for volume 1 by "unrotating" the new voxel locations by original orientation
temp = oldRot_1';
[out.newInterpX_1, out.newInterpY_1, out.newInterpZ_1] = ...
    getRotatedInterpolationGrid(out.locX, out.locY, out.locZ, out.IOP, temp);

%get interpolation points for volume 2 by "unrotating" the new voxel locations by original orientation
temp = oldRot_2';
[out.newInterpX_2, out.newInterpY_2, out.newInterpZ_2] = ...
    getRotatedInterpolationGrid(out.locX, out.locY, out.locZ, out.IOP, temp);

out.voxDim1 = voxDim1;
out.voxDim2 = voxDim2;
out.voxDim3 = voxDim3;

%output final ImagePositionPatient
out.IPP = [squeeze(out.locX(1,1,1:end))'; ...
    squeeze(out.locY(1,1,1:end))'; ...
    squeeze(out.locZ(1,1,1:end))'];

end