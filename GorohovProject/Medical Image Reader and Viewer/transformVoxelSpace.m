function out = transformVoxelSpace(inStruct1,inStruct2)

%this function transforms input volume 1 to voxel space of input volume 2

%check if input volumes have different dimensions
if (size(inStruct1.volumes,3)==1 && size(inStruct2.volumes,3)>1) || ...
        (size(inStruct1.volumes,3)>1 && size(inStruct2.volumes,3)==1)
    %order to transform 3D volume to 2D slice
    if size(inStruct1.volumes,3)==1
        disp('Swapping input order to transform 3D to 2D space')
        temp = inStruct1;
        inStruct1 = inStruct2;
        inStruct2 = temp;
        clear temp
    end
end
%2D data need additional checks
if (size(inStruct1.volumes,3)==1 && size(inStruct2.volumes,3)==1)
    if~all(inStruct1.IOP==inStruct2.IOP)
        disp('Input slices have different orientations and cannot be matched')
        out = [];
        return
    end
end

if ~isfield(inStruct1,'locX')
    [inStruct1.locX,inStruct1.locX,inStruct1.locX] = getVoxelLocs(inStruct1);
end
if ~isfield(inStruct2,'locX')
    [inStruct2.locX,inStruct2.locX,inStruct2.locX] = getVoxelLocs(inStruct2);
end

inStruct1.imSz1 = size(inStruct1.locX,1);
inStruct1.imSz2 = size(inStruct1.locX,2);
inStruct1.imSz3 = size(inStruct1.locX,3);
inStruct2.imSz1 = size(inStruct2.locX,1);
inStruct2.imSz2 = size(inStruct2.locX,2);
inStruct2.imSz3 = size(inStruct2.locX,3);

%start with original structure
out = inStruct1;

%check if volumes already have exactly the same sampling
if all(inStruct1.IOP==inStruct2.IOP) && ...
        (inStruct1.imSz1==inStruct2.imSz1 && inStruct1.imSz2==inStruct2.imSz2 && inStruct1.imSz3==inStruct2.imSz3) && ...
        (all(inStruct1.bedRange1==inStruct2.bedRange1) && all(inStruct1.bedRange2==inStruct2.bedRange2) && all(inStruct1.bedRange3==inStruct2.bedRange3))
    return
end

%now update with new values
out.locX = inStruct2.locX;
out.locY = inStruct2.locY;
out.locZ = inStruct2.locZ;
out.bedRange1 = inStruct2.bedRange1;
out.bedRange2 = inStruct2.bedRange2;
out.bedRange3 = inStruct2.bedRange3;
out.IOP = inStruct2.IOP;
out.IPP = inStruct2.IPP;

rot_1 = [inStruct1.IOP(1:3) inStruct1.IOP(4:6) cross(inStruct1.IOP(1:3),inStruct1.IOP(4:6))];
rot_2 = [inStruct2.IOP(1:3) inStruct2.IOP(4:6) cross(inStruct2.IOP(1:3),inStruct2.IOP(4:6))];

%volume 1
if isfield(inStruct1,'volume')
    vol = inStruct1.volume;
    inStruct1 = rmfield(inStruct1,'volume');
elseif isfield(inStruct1,'volumes')
    vol = inStruct1.volumes;
end
bak = vol(1); %background pixel value
if ~isfield(inStruct1,'locX')
    [inStruct1.locX,inStruct1.locY,inStruct1.locZ] = getVoxelLocs(inStruct1);
end
%volume 2
if ~isfield(inStruct2,'locX')
    [inStruct2.locX,inStruct2.locY,inStruct2.locZ] = getVoxelLocs(inStruct2);
end

%precheck new interpolation points
temp = rotateInterpolationPoints(inStruct1,inStruct2.IOP);

%check for special case: same output voxel space, e.g. same volumes but different orientations
if all(size(temp.newLocX)==size(inStruct2.locX)) &&...
        temp.newLocX(1,1,1)==inStruct2.locX(1,1,1) && ...
        temp.newLocX(end,1,1)==inStruct2.locX(end,1,1) && ...
        temp.newLocX(1,end,1)==inStruct2.locX(1,end,1) && ...
        temp.newLocX(1,1,end)==inStruct2.locX(1,1,end) && ...
        temp.newLocY(1,1,1)==inStruct2.locY(1,1,1) && ...
        temp.newLocY(end,1,1)==inStruct2.locY(end,1,1) && ...
        temp.newLocY(1,end,1)==inStruct2.locY(1,end,1) && ...
        temp.newLocY(1,1,end)==inStruct2.locY(1,1,end) && ...
        temp.newLocZ(1,1,1)==inStruct2.locZ(1,1,1) && ...
        temp.newLocZ(end,1,1)==inStruct2.locZ(end,1,1) && ...
        temp.newLocZ(1,end,1)==inStruct2.locZ(1,end,1) && ...
        temp.newLocZ(1,1,end)==inStruct2.locZ(1,1,end)
    %check for 90 degree rotations, i.e. only axes swapping - no interpolation needed
    prec = 1e9;
    vectTest = round((rot_1'*rot_2)*prec)/prec;
    if all(vectTest(:)==0 | abs(vectTest(:))==1)
        %the new interpolation points are calculated by simply flipping the old interpolation points
        %not by rotating the acutal new voxel locations as is done below (for more complicated rotations)
        if ~isfield(inStruct1,'voxDim1')
            inStruct1.voxDim1 = diff(inStruct1.bedRange1)/size(inStruct1.locX,1);
            inStruct1.voxDim2 = diff(inStruct1.bedRange2)/size(inStruct1.locX,2);
            inStruct1.voxDim3 = diff(inStruct1.bedRange3)/size(inStruct1.locX,3);
        end
        if ~isfield(inStruct2,'voxDim1')
            inStruct2.voxDim1 = diff(inStruct2.bedRange1)/size(inStruct2.locX,1);
            inStruct2.voxDim2 = diff(inStruct2.bedRange2)/size(inStruct2.locX,2);
            inStruct2.voxDim3 = diff(inStruct2.bedRange3)/size(inStruct2.locX,3);
        end
        [interpOutX, interpOutY, interpOutZ] = locFlip(temp.oldInterpX,temp.oldInterpY,temp.oldInterpZ,vectTest,...
            [inStruct1.voxDim1 inStruct1.voxDim2 inStruct1.voxDim3],[inStruct2.voxDim1 inStruct2.voxDim2 inStruct2.voxDim3]);
    else
        %get new points for new space, must interpolate
        [interpOutX, interpOutY, interpOutZ] = ...
            getRotatedInterpolationGrid(inStruct2.locX, inStruct2.locY, inStruct2.locZ, inStruct2.IOP, rot_1');
    end
else
    %get interpolation points for volume 1 by "unrotating" the new voxel locations by original orientation
    [interpOutX, interpOutY, interpOutZ] = ...
        getRotatedInterpolationGrid(inStruct2.locX, inStruct2.locY, inStruct2.locZ, inStruct2.IOP, rot_1');
end

%interpolate
if ((size(inStruct1.volumes,3)==1 && size(inStruct2.volumes,3)>1) || ...
        (size(inStruct1.volumes,3)>1 && size(inStruct2.volumes,3)==1)) || ...
        ~(all(size(temp.oldInterpX)==size(interpOutX)) && ...
        (temp.oldInterpX(1,1,1)==interpOutX(1,1,1) && temp.oldInterpX(end,1,1)==interpOutX(end,1,1) && ...
        temp.oldInterpX(1,end,1)==interpOutX(1,end,1) && temp.oldInterpX(1,1,end)==interpOutX(1,1,end) && ...
        temp.oldInterpY(1,1,1)==interpOutY(1,1,1) && temp.oldInterpY(end,1,1)==interpOutY(end,1,1) && ...
        temp.oldInterpY(1,end,1)==interpOutY(1,end,1) && temp.oldInterpY(1,1,end)==interpOutY(1,1,end) && ...
        temp.oldInterpZ(1,1,1)==interpOutZ(1,1,1) && temp.oldInterpZ(end,1,1)==interpOutZ(end,1,1) && ...
        temp.oldInterpZ(1,end,1)==interpOutZ(1,end,1) && temp.oldInterpZ(1,1,end)==interpOutZ(1,1,end)))
    %data must be single or double
    cst = 0;
    if ~(isa(vol,'single') || isa(vol,'double'))
        origClass = class(vol);
        vol = single(vol);
        cst = 1;
    end
    outVol = interpolateVolume(temp.oldInterpX,temp.oldInterpY,temp.oldInterpZ,vol,interpOutX,interpOutY,interpOutZ);
    %set interpolation NaN values to background
    outVol(isnan(outVol)) = bak;
    if cst
        if islogical(vol)
            outVol = round(outVol);
        end
        eval(['outVol = ' origClass '(outVol);'])
    end
else
    outVol = vol;
end

%set output fields to new values
out.volumes = outVol;
out.imSz1 = size(outVol,1);
out.imSz2 = size(outVol,2);
out.imSz3 = size(outVol,3);

end