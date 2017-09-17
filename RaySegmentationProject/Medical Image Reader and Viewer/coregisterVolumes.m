function out = coregisterVolumes(in1, in2, varargin)

%this function fuses two volumes, making them the same size and dimensions
%inputs are 2 structures from 'readImages' function
%varargin allows option for volume trimming

%throw error if one set is 2D and the other is 3D
if (size(in1.volumes,3)==1 && size(in2.volumes,3)>1) || ...
        (size(in1.volumes,3)>1 && size(in2.volumes,3)==1)
    disp('Inputs have different dimensions and should not be coregistered')
    out = [];
    return
end
%throw error if both sets are 2D and have different orientations
if (size(in1.volumes,3)==1 && size(in2.volumes,3)==1) && ~all(in1.IOP==in2.IOP)
    disp('Input slices have different orientations and cannot be coregistered')
    out = [];
    return
end
%throw warning if data sets are from different scanners
if isfield(in1,'gantryModel') && isfield(in2,'gantryModel')
    if strcmpi(in1.gantryModel(end-2:end),'mCT') || ...
            strcmpi(in1.gantryModel,'1104')
        if ~strcmpi(in2.gantryModel(end-2:end),'mCT') && ...
                ~strcmpi(in2.gantryModel,'1104')
            disp('Gantry model names are different in both data sets')
        end
    elseif strcmpi(in1.gantryModel(end-2:end),'mMR') || ...
            strcmpi(in1.gantryModel,'2008')
        if ~strcmpi(in2.gantryModel(end-2:end),'mMR') && ...
                ~strcmpi(in2.gantryModel,'2008')
            disp('Gantry model names are different in both data sets')
        end
    end
end

h = waitbar(0,'Coregistering volumes...');

%% get matrix of new sampling points

temp = coregisterLocs(in1,in2,varargin{:},'coreg');
clear in1 in2

if exist('h','var') && ishandle(h) && strcmp(get(h,'type'),'figure')
    waitbar(0.5,h,'Interpolating...');
end

%% interpolate over new points

out = coregisterInterp(temp);

if exist('h','var') && ishandle(h) && strcmp(get(h,'type'),'figure')
    waitbar(1,h,'Coregistration done');
end

if exist('h','var') && ishandle(h) && strcmp(get(h,'type'),'figure')
%     pause(0.5)
    close(h)
    drawnow
end

end