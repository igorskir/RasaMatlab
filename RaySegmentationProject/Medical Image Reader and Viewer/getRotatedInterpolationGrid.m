function [outX, outY, outZ] = getRotatedInterpolationGrid(inX, inY, inZ, oldIOP, newRot)

%this function returns a rotated grid of coordinates to be used for all interpolation operations

oldIOP = oldIOP(:); %this is the current orientation
%rotation matrix
oldRot = [oldIOP(1:3) oldIOP(4:6) cross(oldIOP(1:3),oldIOP(4:6))];

%newRot is the relative rotation, NOT the final space

prec = 1e9; %rotation precision
totRot = round(newRot * oldRot * prec) / prec; %this is the new orientation SPACE
oldRot = round(oldRot * prec) / prec;
newRot = round(newRot * prec) / prec;

%check first for simple rotations, not needing point interpolation
if all(newRot(:)==0 | abs(newRot(:))==1)
    [outX, outY, outZ] = simpleRotate(inX, inY, inZ, newRot);
    return
end

imSz1 = size(inX,1);
imSz2 = size(inX,2);
imSz3 = size(inX,3);

%get relative voxel dimensions
temp = oldRot'*[inX(1,1,1) inY(1,1,1) inZ(1,1,1)]';
tempX = oldRot'*[inX(end,1,1) inY(end,1,1) inZ(end,1,1)]';
tempY = oldRot'*[inX(1,end,1) inY(1,end,1) inZ(1,end,1)]';
voxDim1 = (tempX(1) - temp(1))/(imSz1 - 1);
voxDim2 = (tempY(2) - temp(2))/(imSz2 - 1);

%now perform rotation on unrotated grid

%which row and column to input into 'getVoxelLocs'
idx1 = 1;
idx2 = 1;
if voxDim1<0, idx1 = imSz1; end
if voxDim2<0, idx2 = imSz2; end

%assure that matrices are monotonic (probably unnecessary)
tempX = squeeze(double(inX(idx1,idx2,1:end)))';
tempY = squeeze(double(inY(idx1,idx2,1:end)))';
tempZ = squeeze(double(inZ(idx1,idx2,1:end)))';

intX = (tempX(end)-tempX(1))/(imSz3-1);
intY = (tempY(end)-tempY(1))/(imSz3-1);
intZ = (tempZ(end)-tempZ(1))/(imSz3-1);
if intX~=0
    tempX = (-intX*(imSz3-1)/2:intX:intX*(imSz3-1)/2) + mean(tempX);
end
if intY~=0
    tempY = (-intY*(imSz3-1)/2:intY:intY*(imSz3-1)/2) + mean(tempY);
end
if intZ~=0
    tempZ = (-intZ*(imSz3-1)/2:intZ:intZ*(imSz3-1)/2) + mean(tempZ);
end

temp = newRot * [tempX; tempY; tempZ];

if all(totRot(:)==[1;0;0;0;1;0;0;0;1])
    %assure aligned grid in X and Y is not affected by rotation precision
    temp(1,:) = temp(1,1);
    temp(2,:) = temp(2,1);
end

%should use 'abs' and 'flip' method below to exactly match interpolation points
[outX, outY, outZ] = getVoxelLocs(totRot(1:6), temp, ...
    imSz1, imSz2, imSz3, abs(voxDim1), abs(voxDim2));

if voxDim1<0
    try
        outX = flip(outX,1);
        outY = flip(outY,1);
        outZ = flip(outZ,1);
    catch
        outX = flipdim(outX,1);
        outY = flipdim(outY,1);
        outZ = flipdim(outZ,1);
    end
end
if voxDim2<0
    try
        outX = flip(outX,2);
        outY = flip(outY,2);
        outZ = flip(outZ,2);
    catch
        outX = flipdim(outX,2);
        outY = flipdim(outY,2);
        outZ = flipdim(outZ,2);
    end
end

end

function [outX, outY, outZ] = simpleRotate(inX, inY, inZ, newRot)

%rotation matrix should contain only 1's and 0's
prec = 1e9;
newRot = round((newRot)*prec)/prec;
if ~(all(newRot(:)==0 | abs(newRot(:))==1))
    %no reason to get here
    error('This matrix is not a flipping rotation')
end

if abs(newRot(1,1))==1
    outX = inX;
    if abs(newRot(2,2))==1
        %same axes (but possible different directions)
        outY = inY;
        outZ = inZ;
    elseif abs(newRot(3,2))==1
        outY = inZ;
        outZ = inY;
    end
elseif abs(newRot(2,1))==1
    outY = inX;
    if abs(newRot(1,2))==1
        outX = inY;
        outZ = inZ;
    elseif abs(newRot(3,2))==1
        outX = inZ;
        outZ = inY;
    end
elseif abs(newRot(3,1))==1
    outZ = inX;
    if abs(newRot(1,2))==1
        outX = inY;
        outY = inZ;
    elseif abs(newRot(2,2))==1
        outX = inZ;
        outY = inY;
    end
end

if sum(newRot(1,:))==-1
    outX = -outX;
end
if sum(newRot(2,:))==-1
    outY = -outY;
end
if sum(newRot(3,:))==-1
    outZ = -outZ;
end

end