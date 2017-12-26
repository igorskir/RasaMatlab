function out = coregisterInterp(out)

%input for this function is output from 'coregisterLocs'
%to efficiently handle memory, modify input instead of creating new variables

bak1 = out.vol_1(1);
bak2 = out.vol_2(1);
cst1 = 0;
cst2 = 0;

if length(size(out.oldInterpX_1))~=length(size(out.newInterpX_1)) || ...
        ~(all(size(out.oldInterpX_1)==size(out.newInterpX_1)) && ...
        (out.oldInterpX_1(1,1,1)==out.newInterpX_1(1,1,1) && ...
        out.oldInterpX_1(end,1,1)==out.newInterpX_1(end,1,1) && ...
        out.oldInterpX_1(1,end,1)==out.newInterpX_1(1,end,1) && ...
        out.oldInterpX_1(1,1,end)==out.newInterpX_1(1,1,end) && ...
        out.oldInterpY_1(1,1,1)==out.newInterpY_1(1,1,1) && ...
        out.oldInterpY_1(end,1,1)==out.newInterpY_1(end,1,1) && ...
        out.oldInterpY_1(1,end,1)==out.newInterpY_1(1,end,1) && ...
        out.oldInterpY_1(1,1,end)==out.newInterpY_1(1,1,end) && ...
        out.oldInterpZ_1(1,1,1)==out.newInterpZ_1(1,1,1) && ...
        out.oldInterpZ_1(end,1,1)==out.newInterpZ_1(end,1,1) && ...
        out.oldInterpZ_1(1,end,1)==out.newInterpZ_1(1,end,1) && ...
        out.oldInterpZ_1(1,1,end)==out.newInterpZ_1(1,1,end)))
    %data must be single or double
    if ~(isa(bak1,'single') || isa(bak1,'double'))
        origClass1 = class(bak1);
        out.vol_1 = single(out.vol_1);
        cst1 = 1;
    end
    out.vol_1 = interpolateVolume(out.oldInterpX_1,out.oldInterpY_1,out.oldInterpZ_1,...
        out.vol_1,out.newInterpX_1,out.newInterpY_1,out.newInterpZ_1);
end
if length(size(out.oldInterpX_2))~=length(size(out.newInterpX_2)) || ...
        ~(all(size(out.oldInterpX_2)==size(out.newInterpX_2)) && ...
        (out.oldInterpX_2(1,1,1)==out.newInterpX_2(1,1,1) && ...
        out.oldInterpX_2(end,1,1)==out.newInterpX_2(end,1,1) && ...
        out.oldInterpX_2(1,end,1)==out.newInterpX_2(1,end,1) && ...
        out.oldInterpX_2(1,1,end)==out.newInterpX_2(1,1,end) && ...
        out.oldInterpY_2(1,1,1)==out.newInterpY_2(1,1,1) && ...
        out.oldInterpY_2(end,1,1)==out.newInterpY_2(end,1,1) && ...
        out.oldInterpY_2(1,end,1)==out.newInterpY_2(1,end,1) && ...
        out.oldInterpY_2(1,1,end)==out.newInterpY_2(1,1,end) && ...
        out.oldInterpZ_2(1,1,1)==out.newInterpZ_2(1,1,1) && ...
        out.oldInterpZ_2(end,1,1)==out.newInterpZ_2(end,1,1) && ...
        out.oldInterpZ_2(1,end,1)==out.newInterpZ_2(1,end,1) && ...
        out.oldInterpZ_2(1,1,end)==out.newInterpZ_2(1,1,end)))
    %data must be single or double
    if ~(isa(bak2,'single') || isa(bak2,'double'))
        origClass2 = class(bak2);
        out.vol_2 = single(out.vol_2);
        cst2 = 1;
    end
    out.vol_2 = interpolateVolume(out.oldInterpX_2,out.oldInterpY_2,out.oldInterpZ_2,...
        out.vol_2,out.newInterpX_2,out.newInterpY_2,out.newInterpZ_2);
end

%now find NaN values from interpolation
nanIdx_1 = isnan(out.vol_1);
nanIdx_2 = isnan(out.vol_2);
%set NaN values to background
out.vol_1(nanIdx_1) = bak1;
out.vol_2(nanIdx_2) = bak2;
%use NaN indeces to find any complete rows, columns, or slices that are NaN in each volume
nanIdx_1 = nanIdx_1(:,:,:,end);
nanIdx_2 = nanIdx_2(:,:,:,end);
nanIdx1_1 = sum(sum(nanIdx_1,2),3)==size(out.vol_1,2)*size(out.vol_1,3);
nanIdx2_1 = sum(sum(nanIdx_1,1),3)==size(out.vol_1,1)*size(out.vol_1,3);
nanIdx3_1 = sum(sum(nanIdx_1,1),2)==size(out.vol_1,1)*size(out.vol_1,2);
clear nanIdx_1
nanIdx1_2 = sum(sum(nanIdx_2,2),3)==size(out.vol_2,2)*size(out.vol_2,3);
nanIdx2_2 = sum(sum(nanIdx_2,1),3)==size(out.vol_2,1)*size(out.vol_2,3);
nanIdx3_2 = sum(sum(nanIdx_2,1),2)==size(out.vol_2,1)*size(out.vol_2,2);
clear nanIdx_2
keepIdx1 = ~(nanIdx1_1&nanIdx1_2);
keepIdx2 = ~(nanIdx2_1&nanIdx2_2);
keepIdx3 = ~(nanIdx3_1&nanIdx3_2);
%only remove edge values
idx = find(keepIdx1==1);
keepIdx1(idx(1):idx(end)) = 1;
idx = find(keepIdx2==1);
keepIdx2(idx(1):idx(end)) = 1;
idx = find(keepIdx3==1);
keepIdx3(idx(1):idx(end)) = 1;
clear nanIdx1_1 nanIdx1_2 nanIdx2_1 nanIdx2_2 nanIdx3_1 nanIdx3_2
%remove NaN values
if sum(keepIdx1)<length(keepIdx1) || sum(keepIdx2)<length(keepIdx2) || sum(keepIdx3)<length(keepIdx3)
%     temp = out.vol_1(keepIdx1,keepIdx2,keepIdx3);
%     out.vol_1 = temp;
%     clear temp
%     temp = out.vol_2(keepIdx1,keepIdx2,keepIdx3);
%     out.vol_2 = temp;
%     clear temp
%     temp = out.locX(keepIdx1,keepIdx2,keepIdx3);
%     out.locX = temp;
%     clear temp
%     temp = out.locY(keepIdx1,keepIdx2,keepIdx3);
%     out.locY = temp;
%     clear temp
%     temp = out.locZ(keepIdx1,keepIdx2,keepIdx3);
%     out.locZ = temp;
%     clear temp
end
clear keepIdx1 keepIdx2 keepIdx3

%recast if needed
if cst1
    if islogical(bak1)
        out.vol_1 = round(out.vol_1);
    end
    eval(['out.vol_1 = ' origClass1 '(out.vol_1);'])
end
if cst2
    if islogical(bak2)
        out.vol_2 = round(out.vol_2);
    end
    eval(['out.vol_2 = ' origClass2 '(out.vol_2);'])
end

%now calculate bed ranges (along patient axes)
cntr1 = round(size(out.locX,1)/2);
cntr2 = round(size(out.locX,2)/2);
cntr3 = round(size(out.locX,3)/2);
vect1 = out.IOP(1:3);
vect2 = out.IOP(4:6);
vect3 = cross(vect1, vect2);
%bed starts and ends along new axes
bedPts1 = vect1' * [out.locX(1,cntr2,cntr3) out.locY(1,cntr2,cntr3) out.locZ(1,cntr2,cntr3); ...
    out.locX(end,cntr2,cntr3) out.locY(end,cntr2,cntr3) out.locZ(end,cntr2,cntr3)]';
bedPts2 = vect2' * [out.locX(cntr1,1,cntr3) out.locY(cntr1,1,cntr3) out.locZ(cntr1,1,cntr3); ...
    out.locX(cntr1,end,cntr3) out.locY(cntr1,end,cntr3) out.locZ(cntr1,end,cntr3)]';
bedPts3 = vect3' * [out.locX(cntr1,cntr2,1) out.locY(cntr1,cntr2,1) out.locZ(cntr1,cntr2,1); ...
    out.locX(cntr1,cntr2,end) out.locY(cntr1,cntr2,end) out.locZ(cntr1,cntr2,end)]';

%define bed ranges, use 'abs' method
if bedPts1(1)<=bedPts1(2)
    out.bedRange1 = [bedPts1(1)-abs(out.voxDim1)/2 bedPts1(2)+abs(out.voxDim1)/2];
else
    out.bedRange1 = [bedPts1(1)+abs(out.voxDim1)/2 bedPts1(2)-abs(out.voxDim1)/2];
end
if bedPts2(1)<=bedPts2(2)
    out.bedRange2 = [bedPts2(1)-abs(out.voxDim2)/2 bedPts2(2)+abs(out.voxDim2)/2];
else
    out.bedRange2 = [bedPts2(1)+abs(out.voxDim2)/2 bedPts2(2)-abs(out.voxDim2)/2];
end
if bedPts3(1)<=bedPts3(2)
    out.bedRange3 = [bedPts3(1)-abs(out.voxDim3)/2 bedPts3(2)+abs(out.voxDim3)/2];
else
    out.bedRange3 = [bedPts3(1)+abs(out.voxDim3)/2 bedPts3(2)-abs(out.voxDim3)/2];
end

%variables to be shared by both volumes
out.IOP = out.IOP;
out.IPP = [squeeze(out.locX(1,1,1:end))'; ...
    squeeze(out.locY(1,1,1:end))'; ...
    squeeze(out.locZ(1,1,1:end))'];

end