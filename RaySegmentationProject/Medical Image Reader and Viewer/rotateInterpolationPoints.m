function out = rotateInterpolationPoints(in, varargin)

%this function returns old and new points for interpolation (usually different matrix sizes)
%input structure from 'readImages function' and new IOP

coregister = false; %coregistration will not require individual output sampling

if ~isempty(varargin)
    for n=1:length(varargin)
        if isnumeric(varargin{n})
            if length(varargin{n}(:))==9 %rotation format
                newIOP = varargin{n}(1:6)';
                %assure normalized vectors
                newIOP(1:3) = newIOP(1:3)/norm(newIOP(1:3));
                newIOP(4:6) = newIOP(4:6)/norm(newIOP(4:6));
                %check if valid IOP vector
                if any(isnan(newIOP)) || (round(newIOP(1:3)'*newIOP(4:6)*1e6)/1e6)~=0
                    error('Invalid patient orientation input vector');
                end
                %new rotation (aligns to newIOP)
                xr = newIOP(1);
                yr = newIOP(2);
                zr = newIOP(3);
                xc = newIOP(4);
                yc = newIOP(5);
                zc = newIOP(6);
                xs = yr*zc-yc*zr;
                ys = zr*xc-zc*xr;
                zs = xr*yc-xc*yr;
                newRot = [xr xc xs; yr yc ys; zr zc zs];
            elseif length(varargin{n}(:))==6 %IOP format
                newIOP = varargin{:}(:);
                %assure normalized vectors
                newIOP(1:3) = newIOP(1:3)/norm(newIOP(1:3));
                newIOP(4:6) = newIOP(4:6)/norm(newIOP(4:6));
                %check if valid IOP vector
                if any(isnan(newIOP)) || (round(newIOP(1:3)'*newIOP(4:6)*1e6)/1e6)~=0
                    error('Invalid patient orientation input vector');
                end
                %new rotation (aligns to newIOP)
                xr = newIOP(1);
                yr = newIOP(2);
                zr = newIOP(3);
                xc = newIOP(4);
                yc = newIOP(5);
                zc = newIOP(6);
                xs = yr*zc-yc*zr;
                ys = zr*xc-zc*xr;
                zs = xr*yc-xc*yr;
                newRot = [xr xc xs; yr yc ys; zr zc zs];
            end
        end
        if ischar(varargin{n}) && strcmpi(varargin{n}(1),'c')
            coregister = true;
        end
    end
end

%if no orientation vector input, align to scanner coordinates
if ~exist('newRot','var')
    newIOP = [1;0;0;0;1;0];
    newRot = eye(3);
end

if ~isfield(in,'locX')
    [oldLocX,oldLocY,oldLocZ] = getVoxelLocs(in);
else
    oldLocX = in.locX;
    oldLocY = in.locY;
    oldLocZ = in.locZ;
end

oldImSz1 = size(oldLocX,1);
oldImSz2 = size(oldLocX,2);
oldImSz3 = size(oldLocX,3);
oldIOP = in.IOP;

%precision to 0.00001 mm
out.oldVoxDim1 = round(diff(double(in.bedRange1))/oldImSz1*1e5)/1e5;
out.oldVoxDim2 = round(diff(double(in.bedRange2))/oldImSz2*1e5)/1e5;
out.oldVoxDim3 = round(diff(double(in.bedRange3))/oldImSz3*1e5)/1e5;
clear in

%old rotation (aligns to scanner)
xr = oldIOP(1);
yr = oldIOP(2);
zr = oldIOP(3);
xc = oldIOP(4);
yc = oldIOP(5);
zc = oldIOP(6);
xs = yr*zc-yc*zr;
ys = zr*xc-zc*xr;
zs = xr*yc-xc*yr;
oldRot = [xr xc xs; yr yc ys; zr zc zs];

oldFOV1 = out.oldVoxDim1 * oldImSz1;
oldFOV2 = out.oldVoxDim2 * oldImSz2;
oldFOV3 = out.oldVoxDim3 * oldImSz3;

%define new span of resliced image
[newFOV1, newFOV2, newFOV3] = rotateValues(oldFOV1,oldFOV2,oldFOV3,oldIOP,newIOP);
%get new image matrix dimensions
[newImSz1, newImSz2, newImSz3] = rotateValues(oldImSz1,oldImSz2,oldImSz3,oldIOP,newIOP);
clear temp

newImSz1 = round(newImSz1);
newImSz2 = round(newImSz2);
newImSz3 = round(newImSz3);
%new voxel sizes, precision to 0.00001 mm
out.newVoxDim1 = round(newFOV1/newImSz1*1e5)/1e5;
out.newVoxDim2 = round(newFOV2/newImSz2*1e5)/1e5;
out.newVoxDim3 = round(newFOV3/newImSz3*1e5)/1e5;

%get aligned input grid format by unrotation with old IOP
[out.oldInterpX,out.oldInterpY,out.oldInterpZ] = ...
    getRotatedInterpolationGrid(oldLocX,oldLocY,oldLocZ,oldIOP,oldRot');
%assure NDGRID format for interpolation input
[out.oldInterpY,out.oldInterpX,out.oldInterpZ] = ...
    meshgrid(squeeze(out.oldInterpY(1,1:end,1)),squeeze(out.oldInterpX(1:end,1,1)),squeeze(out.oldInterpZ(1,1,1:end)));

%new voxel dimensions
if out.oldInterpX(1,1,1)>out.oldInterpX(2,1,1)
    out.newVoxDim1 = -out.newVoxDim1;
end
if out.oldInterpY(1,1,1)>out.oldInterpY(1,2,1)
    out.newVoxDim2 = -out.newVoxDim2;
end
if out.oldInterpZ(1,1,1)>out.oldInterpZ(1,1,2)
    out.newVoxDim3 = -out.newVoxDim3;
end

%check for 90 degree rotations, i.e. only axes swapping - no interpolation needed
prec = 1e9;
vectTest = round((oldRot'*newRot)*prec)/prec;
if all(vectTest(:)==0 | abs(vectTest(:))==1)
    %check
    if ~all(ismember([newImSz1 newImSz2 newImSz3],[oldImSz1 oldImSz2 oldImSz3]))
        error('Something weird happened')
    end
    %flip the original voxel locations to preserve exact precision
    [out.newLocX, out.newLocY, out.newLocZ] = locFlip(oldLocX,oldLocY,oldLocZ,vectTest,...
        [out.oldVoxDim1 out.oldVoxDim2 out.oldVoxDim3],[out.newVoxDim1 out.newVoxDim2 out.newVoxDim3]);
    %output new rotated interpolation points if not called by coregistering function
    if ~coregister
        %notice that the new interpolation points are calculated by simply flipping the old interpolation points
        %not by rotating the acutal new voxel locations as is done below (for more complicated rotations)
        [out.newInterpX, out.newInterpY, out.newInterpZ] = locFlip(out.oldInterpX,out.oldInterpY,out.oldInterpZ,vectTest,...
            [out.oldVoxDim1 out.oldVoxDim2 out.oldVoxDim3],[out.newVoxDim1 out.newVoxDim2 out.newVoxDim3]);
    end
    return
end

%offset to reshift points after rotation
%round to deal with precision errors when calculating isocenter from different orientations
isocntr = round([mean([oldLocX(1) oldLocX(end)]) ...
    mean([oldLocY(1) oldLocY(end)]) ...
    mean([oldLocZ(1) oldLocZ(end)])]*1e2)/1e2; %precision to 0.01 mm
clear oldLocX oldLocY oldLocZ

tempX = single(-out.newVoxDim1*(newImSz1-1)/2:out.newVoxDim1:out.newVoxDim1*(newImSz1-1)/2);
tempY = single(-out.newVoxDim2*(newImSz2-1)/2:out.newVoxDim2:out.newVoxDim2*(newImSz2-1)/2);
tempZ = single(-out.newVoxDim3*(newImSz3-1)/2:out.newVoxDim3:out.newVoxDim3*(newImSz3-1)/2);

[newLocY,newLocX,newLocZ] = meshgrid(tempY,tempX,tempZ);
clear tempX tempY tempZ

%these are the REAL locations, rotated back to output reference frame
[out.newLocX, out.newLocY, out.newLocZ] = ...
    getRotatedInterpolationGrid(newLocX, newLocY, newLocZ, [1;0;0;0;1;0], newRot);
%this assures proper centering
out.newLocX = out.newLocX + isocntr(1);
out.newLocY = out.newLocY + isocntr(2);
out.newLocZ = out.newLocZ + isocntr(3);

%output new rotated interpolation points if not called by coregistering function
if ~coregister
    %get new points by "unrotating" the new voxel locations by original orientation
    [out.newInterpX, out.newInterpY, out.newInterpZ] = ...
        getRotatedInterpolationGrid(out.newLocX, out.newLocY, out.newLocZ, newIOP, oldRot');
end

%output new IPP
out.newIPP = [squeeze(out.newLocX(1,1,1:end))'; ...
    squeeze(out.newLocY(1,1,1:end))'; ...
    squeeze(out.newLocZ(1,1,1:end))'];

end