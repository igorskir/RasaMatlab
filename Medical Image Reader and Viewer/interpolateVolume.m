function outVol = interpolateVolume(inX,inY,inZ,inVol,outX,outY,outZ)

warning('off','MATLAB:griddedInterpolant:MeshgridEval3DWarnId')

if inX(1)>inX(end)
    try
        inX = flip(inX,1);
        inVol = flip(inVol,1);
    catch
        inX = flipdim(inX,1);
        inVol = flipdim(inVol,1);
    end
end
if inY(1)>inY(end)
    try
        inY = flip(inY,2);
        inVol = flip(inVol,2);
    catch
        inY = flipdim(inY,2);
        inVol = flipdim(inVol,2);
    end
end
if inZ(1)>inZ(end)
    try
        inZ = flip(inZ,3);
        inVol = flip(inVol,3);
    catch
        inZ = flipdim(inZ,3);
        inVol = flipdim(inVol,3);
    end
end

%interpolate
outVol = zeros(size(outX,1),size(outX,2),size(outX,3),size(inVol,4),'like',inVol);
if size(inVol,4)==1
        if size(inVol,3)==1 %2D
            F = griddedInterpolant(inX,inY,inVol,'linear','none');
            outVol = F(outX,outY);
        else %3D
            F = griddedInterpolant(inX,inY,inZ,inVol,'linear','none');
            outVol = F(outX,outY,outZ);
        end
else
    for n=1:size(inVol,4)
        if size(inVol,3)==1 %2D
            F = griddedInterpolant(inX,inY,inVol(:,:,:,n),'linear','none');
            outVol(:,:,:,n) = F(outX,outY);
        else %3D
            F = griddedInterpolant(inX,inY,inZ,inVol(:,:,:,n),'linear','none');
            outVol(:,:,:,n) = F(outX,outY,outZ);
        end
    end
end

end