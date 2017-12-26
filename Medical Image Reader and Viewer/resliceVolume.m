function out = resliceVolume(in, varargin)

%this function reslices volume along transaxial plane defined by IOP
%IOP is vector of cosines formatted like dicom field ImageOrientationPatient
%input is structure from 'readImages' function

showSlices = 0; %hardcoded to not show slicing

%if not orientation vector input, assume to align to scanner coordinates
if ~isempty(varargin) && ...
        (isnumeric(varargin{:}) && ...
        (all(size(varargin{:})==[6 1]) || all(size(varargin{:})==[1 6])))
    out.IOP = varargin{:}(:);
    %assure normalized vectors
    out.IOP(1:3) = out.IOP(1:3)/norm(out.IOP(1:3));
    out.IOP(4:6) = out.IOP(4:6)/norm(out.IOP(4:6));
    %check if valid IOP vector
    if any(isnan(out.IOP)) || (round(out.IOP(1:3)'*out.IOP(4:6)*1e6)/1e6)~=0
        error('Invalid patient orientation input vector');
    end
else %scanner coordinates
    out.IOP = [1;0;0;0;1;0];
end

%allow for dynamic frame input
if isfield(in,'volume')
    oldVolume = in.volume;
elseif isfield(in,'volumes')
    oldVolume = in.volumes;
end
bak = oldVolume(1);

% if size(oldVolume,3)<2
%     error('Cannot reslice single planes');
% end

in = rotateInterpolationPoints(in,out.IOP);

%input data must be single or double
cst = false;
if ~(isa(bak,'single') || isa(bak,'double'))
    origClass = class(oldVolume);
    oldVolume = single(oldVolume);
    cst = true;
end

%now reslice
if showSlices %show slicing
    nFrames = size(oldVolume,4);
    %process all frames
    newImSz3 = size(in.newLocX,3);
    out.volume = zeros(size(in.newLocX,1),size(in.newLocX,2),newImSz3,nFrames,'single');
    f1 = figure('Color','w');
    for N=1:nFrames
        for n=1:newImSz3
            sl = slice(in.oldInterpY,in.oldInterpX,in.oldInterpZ,oldVolume(:,:,:,N),...
                in.newInterpY(:,:,n),in.newInterpX(:,:,n),in.newInterpZ(:,:,n));
            set(sl,'EdgeColor','none')
            out.volume(:,:,n,N) = get(sl,'CData');
            axis equal
            view(-44,30)
            xlim([in.oldInterpY(1,1,1) in.oldInterpY(1,end,1)])
            ylim([in.oldInterpX(1,1,1) in.oldInterpX(end,1,1)])
            zlim([in.oldInterpZ(1,1,1) in.oldInterpZ(1,1,end)])
            drawnow
        end
    end
    close(f1)
else
    out.volume = interpolateVolume(in.oldInterpX,in.oldInterpY,in.oldInterpZ,...
        oldVolume,in.newInterpX,in.newInterpY,in.newInterpZ);
end

%now find any NaN values from interpolation
nanIdx = isnan(out.volume);
%set NaN values to background
out.volume(nanIdx) = bak;
%use NaN indeces to find any complete rows, columns, or slices that are entirely NaN
nanIdx = nanIdx(:,:,:,end);
nanIdx1 = sum(sum(nanIdx,2),3)==size(out.volume,2)*size(out.volume,3);
nanIdx2 = sum(sum(nanIdx,1),3)==size(out.volume,1)*size(out.volume,3);
nanIdx3 = sum(sum(nanIdx,1),2)==size(out.volume,1)*size(out.volume,2);
%only remove edge values
idx = find(nanIdx1==0);
if ~isempty(idx)
    nanIdx1(idx(1):idx(end)) = 0;
end
idx = find(nanIdx2==0);
if ~isempty(idx)
    nanIdx2(idx(1):idx(end)) = 0;
end
idx = find(nanIdx3==0);
if ~isempty(idx)
    nanIdx3(idx(1):idx(end)) = 0;
end
clear nanIdx
%remove NaN values
if sum(nanIdx1)>0 || sum(nanIdx2)>0 || sum(nanIdx3)>0
%     if mod(sum(nanIdx1),2)==1 || mod(sum(nanIdx2),2)==1 || mod(sum(nanIdx3),2)==1
%         warning('Interpolation matrices are not centered')
%     end
%     %exact coregistration of identical data can be affected if the output matrix size is changed
%     temp = out.volume(~(nanIdx1),~(nanIdx2),~(nanIdx3),ones(size(out.volume,4),1));
%     out.volume = temp;
%     clear temp
%     temp = in.newLocX(~(nanIdx1),~(nanIdx2),~(nanIdx3));
%     in.newLocX = temp;
%     clear temp
%     temp = in.newLocY(~(nanIdx1),~(nanIdx2),~(nanIdx3));
%     in.newLocY = temp;
%     clear temp
%     temp = in.newLocZ(~(nanIdx1),~(nanIdx2),~(nanIdx3));
%     in.newLocZ = temp;
%     clear temp
end
clear nanIdx1 nanIdx2 nanIdx3

%recast if needed
if cst
    if islogical(bak)
        out.volume = round(out.volume);
    end
    eval(['out.volume = ' origClass '(out.volume);'])
end

out.locX = in.newLocX;
out.locY = in.newLocY;
out.locZ = in.newLocZ;
out.IPP = [squeeze(out.locX(1,1,1:end))'; ...
    squeeze(out.locY(1,1,1:end))'; ...
    squeeze(out.locZ(1,1,1:end))'];

%adjust image sizes
newImSz1 = size(in.newLocX,1);
newImSz2 = size(in.newLocX,2);
newImSz3 = size(in.newLocX,3);

%get bed ranges
% temp1 = [in.newLocX(1,1,1) in.newLocY(1,1,1) in.newLocZ(1,1,1)];
% temp2 = [in.newLocX(end,1,1) in.newLocY(end,1,1) in.newLocZ(end,1,1)];
% out.voxDim1 = norm(temp1-temp2)/(newImSz1-1); %column spacing
% temp2 = [in.newLocX(1,end,1) in.newLocY(1,end,1) in.newLocZ(1,end,1)];
% out.voxDim2 = norm(temp1-temp2)/(newImSz2-1); %row spacing
% temp2 = [in.newLocX(1,1,end) in.newLocY(1,1,end) in.newLocZ(1,1,end)];
% out.voxDim3 = norm(temp1-temp2)/(newImSz3-1); %slice spacing
out.voxDim1 = in.newVoxDim1;
out.voxDim2 = in.newVoxDim2;
out.voxDim3 = in.newVoxDim3;

%calculate bed ranges (along patient axes)
cntr1 = round(newImSz1/2);
cntr2 = round(newImSz2/2);
cntr3 = round(newImSz3/2);
vect1 = out.IOP(1:3);
vect2 = out.IOP(4:6);
vect3 = cross(vect1,vect2);
%bed starts and ends along new axes
bedPts1 = vect1'*[in.newLocX(1,cntr2,cntr3) in.newLocY(1,cntr2,cntr3) in.newLocZ(1,cntr2,cntr3); ...
    in.newLocX(end,cntr2,cntr3) in.newLocY(end,cntr2,cntr3) in.newLocZ(end,cntr2,cntr3)]';
bedPts2 = vect2'*[in.newLocX(cntr1,1,cntr3) in.newLocY(cntr1,1,cntr3) in.newLocZ(cntr1,1,cntr3); ...
    in.newLocX(cntr1,end,cntr3) in.newLocY(cntr1,end,cntr3) in.newLocZ(cntr1,end,cntr3)]';
bedPts3 = vect3'*[in.newLocX(cntr1,cntr2,1) in.newLocY(cntr1,cntr2,1) in.newLocZ(cntr1,cntr2,1); ...
    in.newLocX(cntr1,cntr2,end) in.newLocY(cntr1,cntr2,end) in.newLocZ(cntr1,cntr2,end)]';
%define bed ranges, must use 'abs' method
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

end