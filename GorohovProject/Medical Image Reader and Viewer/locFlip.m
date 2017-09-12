function [outX, outY, outZ] = locFlip(inX,inY,inZ,rot,oldVoxDim,newVoxDim)

%this function handles simple rotations (e.g. 90 degrees) by permutation and flipping

%rotation matrix should contain only 1's and 0's
prec = 1e9;
rot = round((rot)*prec)/prec;
if ~(all(rot(:)==0 | abs(rot(:))==1))
    error('This matrix is not a flipping rotation')
end

permIdx = [find(abs(rot(:,1))==1) find(abs(rot(:,2))==1) find(abs(rot(:,3))==1)];
flipIdx = [find(abs(rot(1,:))==1) find(abs(rot(2,:))==1) find(abs(rot(3,:))==1)];

outX = permute(inX,permIdx);
outY = permute(inY,permIdx);
outZ = permute(inZ,permIdx);

try
    %rotational axes reversal
    if sum(rot(1,:))==-1
        outX = flip(outX,flipIdx(1));
        outY = flip(outY,flipIdx(1));
        outZ = flip(outZ,flipIdx(1));
    end
    if sum(rot(2,:))==-1
        outX = flip(outX,flipIdx(2));
        outY = flip(outY,flipIdx(2));
        outZ = flip(outZ,flipIdx(2));
    end
    if sum(rot(3,:))==-1
        outX = flip(outX,flipIdx(3));
        outY = flip(outY,flipIdx(3));
        outZ = flip(outZ,flipIdx(3));
    end
    %account for relative axes directions
    if newVoxDim(1)/norm(newVoxDim(1))~=oldVoxDim(permIdx(1))/norm(oldVoxDim(permIdx(1)))
        outX = flip(outX,1);
        outY = flip(outY,1);
        outZ = flip(outZ,1);
    end
    if newVoxDim(2)/norm(newVoxDim(2))~=oldVoxDim(permIdx(2))/norm(oldVoxDim(permIdx(2)))
        outX = flip(outX,2);
        outY = flip(outY,2);
        outZ = flip(outZ,2);
    end
    if newVoxDim(3)/norm(newVoxDim(3))~=oldVoxDim(permIdx(3))/norm(oldVoxDim(permIdx(3)))
        outX = flip(outX,3);
        outY = flip(outY,3);
        outZ = flip(outZ,3);
    end
catch
    if sum(rot(1,:))==-1
        outX = flipdim(outX,flipIdx(1));
        outY = flipdim(outY,flipIdx(1));
        outZ = flipdim(outZ,flipIdx(1));
    end
    if sum(rot(2,:))==-1
        outX = flipdim(outX,flipIdx(2));
        outY = flipdim(outY,flipIdx(2));
        outZ = flipdim(outZ,flipIdx(2));
    end
    if sum(rot(3,:))==-1
        outX = flipdim(outX,flipIdx(3));
        outY = flipdim(outY,flipIdx(3));
        outZ = flipdim(outZ,flipIdx(3));
    end
    %account for relative axes directions
    if newVoxDim(1)/norm(newVoxDim(1))~=oldVoxDim(permIdx(1))/norm(oldVoxDim(permIdx(1)))
        outX = flipdim(outX,1);
        outY = flipdim(outY,1);
        outZ = flipdim(outZ,1);
    end
    if newVoxDim(2)/norm(newVoxDim(2))~=oldVoxDim(permIdx(2))/norm(oldVoxDim(permIdx(2)))
        outX = flipdim(outX,2);
        outY = flipdim(outY,2);
        outZ = flipdim(outZ,2);
    end
    if newVoxDim(3)/norm(newVoxDim(3))~=oldVoxDim(permIdx(3))/norm(oldVoxDim(permIdx(3)))
        outX = flipdim(outX,3);
        outY = flipdim(outY,3);
        outZ = flipdim(outZ,3);
    end
end

end