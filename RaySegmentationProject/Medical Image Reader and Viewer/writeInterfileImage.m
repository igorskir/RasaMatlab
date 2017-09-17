function writeInterfileImage(varargin)

%inputs are structure from 'readImages' function and image filename

warning('off','all')

%check core architecture for file path character
test = computer;
if strcmpi(test(1:2),'PC')
    pathChar = '\';
else
    pathChar = '/';
end

[S, imageFile] = parseInputs(varargin{:});

if ~isfield(S,'volume') && ~isfield(S,'volumes')
    error('Input structure should contain a single volume set')
end
%only allow nuclear medicine data to be written to interfile
% if (isfield(S,'modality') && strcmpi(S.modality,'MR')) || ...
%         (isfield(S,'modality') && strcmpi(S.modality,'CT'))
%     error('Interfile format is only suitable for nuclear medicine data.')
% end
if isfield(S,'IOP')
    if ~all((abs(S.IOP)==1)+(abs(S.IOP)==0))
        error('Interfile format is not suitable for rotated volumes')
    end
    IOPcross = cross(S.IOP(1:3),S.IOP(4:6));
    if abs(IOPcross(3))~=1
        error('Unknown patient orientation')
    end
end
if ~isfield(S,'patientOrientation')
    error('Patient orientation is not consistent if volume has been manually rotated')
end

if ~strcmpi(imageFile(end-1:end),'.v'), imageFile = [imageFile '.v']; end
%get folder to make sure it exists
imageDir = regexp(imageFile,pathChar,'split');
try
    imageDir = strjoin(imageDir(1:end-1),pathChar);
catch
    temp = imageDir;
    for n=1:numel(temp)-1
        if n==1
            imageDir = temp{n};
        else
            imageDir = [imageDir pathChar temp{n}];
        end
    end
end

mkdir(imageDir)

%check for previous files
prevFiles = dir(imageDir);
prevFiles(strcmpi({prevFiles.name},'.')) = [];
prevFiles(strcmpi({prevFiles.name},'..')) = [];
% %check for previous 'BAK' folder
% idx = strcmp({prevFiles.name},'BAK');
% if ~isempty(idx)
%     if(prevFiles(idx).isdir)
%         prevFiles(idx) = [];
%     end
% end

%only files
prevFiles([prevFiles.isdir]) = [];
if ~isempty(prevFiles)
    warning('Previous files found in image directory - moving them to folder ''BAK''')
    mkdir([imageDir pathChar 'BAK'])
    for n=1:length(prevFiles)
        movefile([imageDir pathChar prevFiles(n).name],[imageDir pathChar 'BAK' pathChar prevFiles(n).name])
    end
end

if isfield(S,'volume')
    nFrames = size(S.volume,4);
elseif isfield(S,'volumes')
    nFrames = size(S.volumes,4);
else
    error('Input structure does not contain volume data')
end

if nFrames==1
    disp('Writing Interfile image...')
    h = waitbar(0,'Writing Interfile image...');
else
    disp(['Writing ' num2str(nFrames) ' Interfile images...'])
    h = waitbar(0,['Writing ' num2str(nFrames) ' Interfile images...']);
end

meta = initiateInterfileStruct;
meta = getStaticInterfileFields(meta, S);

%correct volume orientation and make consistent with interfile convention
%the following is somewhat cumbersome because corrections are different between mCT and mMR
if isfield(S,'gantryModel') && ...
        (strcmpi(S.gantryModel(end-2:end),'mCT') || strcmpi(S.gantryModel,'1104'))
    if isfield(S,'patientOrientation')
        %mCT interfile IOP may not reflect actual patientOrientation
        switch S.patientOrientation
            case 'HFS'
                if all(S.IOP==[1;0;0;0;1;0])
                    %do nothing
                else
                    if exist('h','var') && ishandle(h) && strcmp(get(h,'type'),'figure')
                        close(h)
                    end
                    error('Patient orientations are inconsistent')
                end
            case 'FFS'
                if all(S.IOP==[1;0;0;0;1;0]) %from interfile and dicom
                    %flip bed start and end and reshift new z offset
                    temp = num2str(-str2double(meta.hBedStart{:}));
                    meta.hBedStart{:} = num2str(-str2double(meta.hBedEnd{:}));
                    meta.hBedEnd{:} = temp;
                else
                    if exist('h','var') && ishandle(h) && strcmp(get(h,'type'),'figure')
                        close(h)
                    end
                    error('Patient orientations are inconsistent')
                end
            case 'HFP'
                if all(S.IOP==[-1;0;0;0;-1;0])
                    %redo vertical bed start
                    meta.vBedStart{:} = num2str(-mean(S.bedRange2)+str2double(meta.yOffset{:}));
                else
                    if exist('h','var') && ishandle(h) && strcmp(get(h,'type'),'figure')
                        close(h)
                    end
                    error('Patient orientations are inconsistent')
                end
            case 'FFP'
                if all(S.IOP==[-1;0;0;0;-1;0])
                    %redo vertical bed start
                    meta.vBedStart{:} = num2str(-mean(S.bedRange2)+str2double(meta.yOffset{:}));
                    %flip bed start and end and reshift new z offset
                    temp = num2str(-str2double(meta.hBedStart{:}));
                    meta.hBedStart{:} = num2str(-str2double(meta.hBedEnd{:}));
                    meta.hBedEnd{:} = temp;
                else
                    if exist('h','var') && ishandle(h) && strcmp(get(h,'type'),'figure')
                        close(h)
                    end
                    error('Patient orientations are inconsistent')
                end
            case 'HFDL'
                if all(S.IOP==[0;-1;0;1;0;0])
                    %redo vertical bed start
                    meta.vBedStart{:} = num2str(-mean(S.bedRange2)+str2double(meta.yOffset{:}));
                else
                    if exist('h','var') && ishandle(h) && strcmp(get(h,'type'),'figure')
                        close(h)
                    end
                    error('Patient orientations are inconsistent')
                end
            case 'FFDL'
                if all(S.IOP==[0;-1;0;1;0;0])
                    %redo vertical bed start
                    meta.vBedStart{:} = num2str(-mean(S.bedRange2)+str2double(meta.yOffset{:}));
                    %flip bed start and end and reshift new z offset
                    temp = num2str(-str2double(meta.hBedStart{:}));
                    meta.hBedStart{:} = num2str(-str2double(meta.hBedEnd{:}));
                    meta.hBedEnd{:} = temp;
                else
                    if exist('h','var') && ishandle(h) && strcmp(get(h,'type'),'figure')
                        close(h)
                    end
                    error('Patient orientations are inconsistent')
                end
            case 'HFDR'
                if all(S.IOP==[0;1;0;-1;0;0])
                    %redo vertical bed start
                    meta.vBedStart{:} = num2str(-mean(S.bedRange2)+str2double(meta.yOffset{:}));
                else
                    if exist('h','var') && ishandle(h) && strcmp(get(h,'type'),'figure')
                        close(h)
                    end
                    error('Patient orientations are inconsistent')
                end
            case 'FFDR'
                if all(S.IOP==[0;1;0;-1;0;0])
                    %redo vertical bed start
                    meta.vBedStart{:} = num2str(-mean(S.bedRange2)+str2double(meta.yOffset{:}));
                    %flip bed start and end and reshift new z offset
                    temp = num2str(-str2double(meta.hBedStart{:}));
                    meta.hBedStart{:} = num2str(-str2double(meta.hBedEnd{:}));
                    meta.hBedEnd{:} = temp;
                else
                    if exist('h','var') && ishandle(h) && strcmp(get(h,'type'),'figure')
                        close(h)
                    end
                    error('Patient orientations are inconsistent')
                end
            otherwise
                if exist('h','var') && ishandle(h) && strcmp(get(h,'type'),'figure')
                    close(h)
                end
                error('Unknown patient orientation')
        end
    else
        %         %currently limited to non-manually rotated volumes
        %         switch num2str(round(S.IOP'*1e3)/1e3)
        %             case '1  0  0  0  1  0'
        %                 S.patientOrientation = 'HFS';
        %             case '-1  0  0  0  1  0'
        %                 S.patientOrientation = 'FFS';
        %             case '-1  0  0  0 -1  0'
        %                 S.patientOrientation = 'HFP';
        %             case '1  0  0  0 -1  0'
        %                 S.patientOrientation = 'FFP';
        %             case '0 -1  0  1  0  0'
        %                 S.patientOrientation = 'HFDL';
        %             case '0  1  0  1  0  0'
        %                 S.patientOrientation = 'FFDL';
        %             case '0  1  0 -1  0  0'
        %                 S.patientOrientation = 'HFDR';
        %             case '0 -1  0 -1  0  0'
        %                 S.patientOrientation = 'FFDR';
        %             otherwise
        %                 error('This shouldn''t have happened')
        %         end
    end
elseif isfield(meta,'gantryModel') && strcmpi(meta.gantryModel{:},'2008')
    if isfield(S,'patientOrientation')
        %mMR interfile always uses correct IOP for each patient orientation
        switch S.patientOrientation
            case 'HFS'
                if all(S.IOP==[1;0;0;0;1;0])
                    %do nothing
                else
                    if exist('h','var') && ishandle(h) && strcmp(get(h,'type'),'figure')
                        close(h)
                    end
                    error('Patient orientations are inconsistent')
                end
            case 'FFS'
                if all(S.IOP==[-1;0;0;0;1;0])
                    %do nothing
                elseif all(S.IOP==[1;0;0;0;1;0]) %if from dicom
                    S.IOP = [-1;0;0;0;1;0];
                    IOPcross = cross(S.IOP(1:3),S.IOP(4:6));
                    if isfield(S,'volume')
                        try
                            S.volume = flip(S.volume,1);
                            S.volume = flip(S.volume,3);
                        catch
                            S.volume = flipdim(S.volume,1);
                            S.volume = flipdim(S.volume,3);
                        end
                    elseif isfield(S,'volumes')
                        try
                            S.volumes = flip(S.volumes,1);
                            S.volumes = flip(S.volumes,3);
                        catch
                            S.volumes = flipdim(S.volumes,1);
                            S.volumes = flipdim(S.volumes,3);
                        end
                    end
                    %unshift old z offset
                    meta.hBedStart{:} = num2str(str2double(meta.hBedStart{:})-str2double(meta.zOffset{:}));
                    meta.hBedEnd{:} = num2str(str2double(meta.hBedEnd{:})-str2double(meta.zOffset{:}));
                    %redo z offset
                    meta.zOffset{:} = num2str(str2double(meta.zOffset{:})*IOPcross(3));
                    %flip bed start and end and reshift new z offset
                    temp = num2str(str2double(meta.hBedStart{:})+str2double(meta.zOffset{:}));
                    meta.hBedStart{:} = num2str(str2double(meta.hBedEnd{:})+str2double(meta.zOffset{:}));
                    meta.hBedEnd{:} = temp;
                else
                    if exist('h','var') && ishandle(h) && strcmp(get(h,'type'),'figure')
                        close(h)
                    end
                    error('Patient orientations are inconsistent')
                end
            case 'HFP'
                if all(S.IOP==[-1;0;0;0;-1;0])
                    %do nothing
                elseif all(S.IOP==[1;0;0;0;1;0]) %if from dicom
                    S.IOP = [-1;0;0;0;-1;0];
                    IOPcross = cross(S.IOP(1:3),S.IOP(4:6));
                    if isfield(S,'volume')
                        try
                            S.volume = flip(S.volume,1);
                            S.volume = flip(S.volume,2);
                        catch
                            S.volume = flipdim(S.volume,1);
                            S.volume = flipdim(S.volume,2);
                        end
                    elseif isfield(S,'volumes')
                        try
                            S.volumes = flip(S.volumes,1);
                            S.volumes = flip(S.volumes,2);
                        catch
                            S.volumes = flipdim(S.volumes,1);
                            S.volumes = flipdim(S.volumes,2);
                        end
                    end
                    %unshift old y offset
                    meta.vBedStart{:} = num2str(str2double(meta.vBedStart{:})-str2double(meta.yOffset{:}));
                    %redo offsets that depend on IOP
                    meta.xOffset{:} = num2str((-1 * (1.02334 + str2double(meta.voxDim2{:})/2) * S.IOP(4)) + ... % HFDL, HFRL, FFDL, or FFDR
                        (-1 * (-0.1450 + (str2double(meta.voxDim1{:})/2 * S.IOP(5))) * S.IOP(1) * S.IOP(5))); % HFS, FFS, HFP, or FFP
                    meta.yOffset{:} = num2str((-1 * (-0.1450 + str2double(meta.voxDim1{:})/2) * S.IOP(2)) + ... % HFDL, HFRL, FFDL, or FFDR
                        (-1 * (1.02334 + (str2double(meta.voxDim2{:})/2 * S.IOP(5))) * S.IOP(5) * S.IOP(5))); % HFS, FFS, HFP, or FFP
                    %reshift new y offset
                    meta.vBedStart{:} = num2str(str2double(meta.vBedStart{:})+str2double(meta.yOffset{:}));
                else
                    if exist('h','var') && ishandle(h) && strcmp(get(h,'type'),'figure')
                        close(h)
                    end
                    error('Patient orientations are inconsistent')
                end
            case 'FFP'
                if all(S.IOP==[1;0;0;0;-1;0])
                    %do nothing
                elseif all(S.IOP==[1;0;0;0;1;0])
                    S.IOP = [1;0;0;0;-1;0];
                    IOPcross = cross(S.IOP(1:3),S.IOP(4:6));
                    if isfield(S,'volume')
                        try
                            S.volume = flip(S.volume,2);
                            S.volume = flip(S.volume,3);
                        catch
                            S.volume = flipdim(S.volume,2);
                            S.volume = flipdim(S.volume,3);
                        end
                    elseif isfield(S,'volumes')
                        try
                            S.volumes = flip(S.volumes,2);
                            S.volumes = flip(S.volumes,3);
                        catch
                            S.volumes = flipdim(S.volumes,2);
                            S.volumes = flipdim(S.volumes,3);
                        end
                    end
                    %unshift old y offset
                    meta.vBedStart{:} = num2str(str2double(meta.vBedStart{:})-str2double(meta.yOffset{:}));
                    %redo offsets that depend on IOP
                    meta.xOffset{:} = num2str((-1 * (1.02334 + str2double(meta.voxDim2{:})/2) * S.IOP(4)) + ... % HFDL, HFRL, FFDL, or FFDR
                        (-1 * (-0.1450 + (str2double(meta.voxDim1{:})/2 * S.IOP(5))) * S.IOP(1) * S.IOP(5))); % HFS, FFS, HFP, or FFP
                    meta.yOffset{:} = num2str((-1 * (-0.1450 + str2double(meta.voxDim1{:})/2) * S.IOP(2)) + ... % HFDL, HFRL, FFDL, or FFDR
                        (-1 * (1.02334 + (str2double(meta.voxDim2{:})/2 * S.IOP(5))) * S.IOP(5) * S.IOP(5))); % HFS, FFS, HFP, or FFP
                    %reshift new y offset
                    meta.vBedStart{:} = num2str(str2double(meta.vBedStart{:})+str2double(meta.yOffset{:}));
                    %unshift old z offset
                    meta.hBedStart{:} = num2str(str2double(meta.hBedStart{:})-str2double(meta.zOffset{:}));
                    meta.hBedEnd{:} = num2str(str2double(meta.hBedEnd{:})-str2double(meta.zOffset{:}));
                    %redo z offset
                    meta.zOffset{:} = num2str(str2double(meta.zOffset{:})*IOPcross(3));
                    %flip bed start and end and reshift new z offset
                    temp = num2str(str2double(meta.hBedStart{:})+str2double(meta.zOffset{:}));
                    meta.hBedStart{:} = num2str(str2double(meta.hBedEnd{:})+str2double(meta.zOffset{:}));
                    meta.hBedEnd{:} = temp;
                else
                    if exist('h','var') && ishandle(h) && strcmp(get(h,'type'),'figure')
                        close(h)
                    end
                    error('Patient orientations are inconsistent')
                end
            case 'HFDL'
                if all(S.IOP==[0;-1;0;1;0;0])
                    %do nothing
                else
                    if exist('h','var') && ishandle(h) && strcmp(get(h,'type'),'figure')
                        close(h)
                    end
                    error('Patient orientations are inconsistent')
                end
            case 'FFDL'
                if all(S.IOP==[0;1;0;1;0;0])
                    %do nothing
                elseif all(S.IOP==[0;-1;0;1;0;0])
                    S.IOP = [0;1;0;1;0;0];
                    IOPcross = cross(S.IOP(1:3),S.IOP(4:6));
                    if isfield(S,'volume')
                        try
                            S.volume = flip(S.volume,1);
                            S.volume = flip(S.volume,3);
                        catch
                            S.volume = flipdim(S.volume,1);
                            S.volume = flipdim(S.volume,3);
                        end
                    elseif isfield(S,'volumes')
                        try
                            S.volumes = flip(S.volumes,1);
                            S.volumes = flip(S.volumes,3);
                        catch
                            S.volumes = flipdim(S.volumes,1);
                            S.volumes = flipdim(S.volumes,3);
                        end
                    end
                    %unshift old y offset
                    meta.vBedStart{:} = num2str(str2double(meta.vBedStart{:})-str2double(meta.yOffset{:}));
                    %redo offsets that depend on IOP
                    meta.xOffset{:} = num2str((-1 * (1.02334 + str2double(meta.voxDim2{:})/2) * S.IOP(4)) + ... % HFDL, HFRL, FFDL, or FFDR
                        (-1 * (-0.1450 + (str2double(meta.voxDim1{:})/2 * S.IOP(5))) * S.IOP(1) * S.IOP(5))); % HFS, FFS, HFP, or FFP
                    meta.yOffset{:} = num2str((-1 * (-0.1450 + str2double(meta.voxDim1{:})/2) * S.IOP(2)) + ... % HFDL, HFRL, FFDL, or FFDR
                        (-1 * (1.02334 + (str2double(meta.voxDim2{:})/2 * S.IOP(5))) * S.IOP(5) * S.IOP(5))); % HFS, FFS, HFP, or FFP
                    %reshift new y offset
                    meta.vBedStart{:} = num2str(str2double(meta.vBedStart{:})+str2double(meta.yOffset{:}));
                    %unshift old z offset
                    meta.hBedStart{:} = num2str(str2double(meta.hBedStart{:})-str2double(meta.zOffset{:}));
                    meta.hBedEnd{:} = num2str(str2double(meta.hBedEnd{:})-str2double(meta.zOffset{:}));
                    %redo z offset
                    meta.zOffset{:} = num2str(str2double(meta.zOffset{:})*IOPcross(3));
                    %flip bed start and end and reshift new z offset
                    temp = num2str(str2double(meta.hBedStart{:})+str2double(meta.zOffset{:}));
                    meta.hBedStart{:} = num2str(str2double(meta.hBedEnd{:})+str2double(meta.zOffset{:}));
                    meta.hBedEnd{:} = temp;
                else
                    if exist('h','var') && ishandle(h) && strcmp(get(h,'type'),'figure')
                        close(h)
                    end
                    error('Patient orientations are inconsistent')
                end
            case 'HFDR'
                if all(S.IOP==[0;1;0;-1;0;0])
                    %do nothing
                else
                    if exist('h','var') && ishandle(h) && strcmp(get(h,'type'),'figure')
                        close(h)
                    end
                    error('Patient orientations are inconsistent')
                end
            case 'FFDR'
                if all(S.IOP==[0;-1;0;-1;0;0])
                    %do nothing
                elseif all(S.IOP==[0;1;0;-1;0;0])
                    S.IOP = [0;-1;0;-1;0;0];
                    IOPcross = cross(S.IOP(1:3),S.IOP(4:6));
                    if isfield(S,'volume')
                        try
                            S.volume = flip(S.volume,1);
                            S.volume = flip(S.volume,3);
                        catch
                            S.volume = flipdim(S.volume,1);
                            S.volume = flipdim(S.volume,3);
                        end
                    elseif isfield(S,'volumes')
                        try
                            S.volumes = flip(S.volumes,1);
                            S.volumes = flip(S.volumes,3);
                        catch
                            S.volumes = flipdim(S.volumes,1);
                            S.volumes = flipdim(S.volumes,3);
                        end
                    end
                    %unshift old y offset
                    meta.vBedStart{:} = num2str(str2double(meta.vBedStart{:})-str2double(meta.yOffset{:}));
                    %redo offsets that depend on IOP
                    meta.xOffset{:} = num2str((-1 * (1.02334 + str2double(meta.voxDim2{:})/2) * S.IOP(4)) + ... % HFDL, HFRL, FFDL, or FFDR
                        (-1 * (-0.1450 + (str2double(meta.voxDim1{:})/2 * S.IOP(5))) * S.IOP(1) * S.IOP(5))); % HFS, FFS, HFP, or FFP
                    meta.yOffset{:} = num2str((-1 * (-0.1450 + str2double(meta.voxDim1{:})/2) * S.IOP(2)) + ... % HFDL, HFRL, FFDL, or FFDR
                        (-1 * (1.02334 + (str2double(meta.voxDim2{:})/2 * S.IOP(5))) * S.IOP(5) * S.IOP(5))); % HFS, FFS, HFP, or FFP
                    %reshift new y offset
                    meta.vBedStart{:} = num2str(str2double(meta.vBedStart{:})+str2double(meta.yOffset{:}));
                    %unshift old z offset
                    meta.hBedStart{:} = num2str(str2double(meta.hBedStart{:})-str2double(meta.zOffset{:}));
                    meta.hBedEnd{:} = num2str(str2double(meta.hBedEnd{:})-str2double(meta.zOffset{:}));
                    %redo z offset
                    meta.zOffset{:} = num2str(str2double(meta.zOffset{:})*IOPcross(3));
                    %flip bed start and end and reshift new z offset
                    temp = num2str(str2double(meta.hBedStart{:})+str2double(meta.zOffset{:}));
                    meta.hBedStart{:} = num2str(str2double(meta.hBedEnd{:})+str2double(meta.zOffset{:}));
                    meta.hBedEnd{:} = temp;
                else
                    if exist('h','var') && ishandle(h) && strcmp(get(h,'type'),'figure')
                        close(h)
                    end
                    error('Patient orientations are inconsistent')
                end
            otherwise
                if exist('h','var') && ishandle(h) && strcmp(get(h,'type'),'figure')
                    close(h)
                end
                error('Unknown patient orientation')
        end
    else
        %         %currently limited to non-manually rotated volumes
        %         switch num2str(round(S.IOP'*1e3)/1e3)
        %             case '1  0  0  0  1  0'
        %                 S.patientOrientation = 'HFS';
        %             case '-1  0  0  0  1  0'
        %                 S.patientOrientation = 'FFS';
        %             case '-1  0  0  0 -1  0'
        %                 S.patientOrientation = 'HFP';
        %             case '1  0  0  0 -1  0'
        %                 S.patientOrientation = 'FFP';
        %             case '0 -1  0  1  0  0'
        %                 S.patientOrientation = 'HFDL';
        %             case '0  1  0  1  0  0'
        %                 S.patientOrientation = 'FFDL';
        %             case '0  1  0 -1  0  0'
        %                 S.patientOrientation = 'HFDR';
        %             case '0 -1  0 -1  0  0'
        %                 S.patientOrientation = 'FFDR';
        %             otherwise
        %                 error('This shouldn''t have happened')
        %         end
    end
    %for mMR
    if str2double(meta.hBedStart{:})*IOPcross(3)<str2double(meta.hBedEnd{:})*IOPcross(3)
        if isfield(S,'volume')
            try
                S.volume = flip(S.volume,3);
            catch
                S.volume = flipdim(S.volume,3);
            end
        elseif isfield(S,'volumes')
            try
                S.volumes = flip(S.volumes,3);
            catch
                S.volumes = flipdim(S.volumes,3);
            end
        end
    end
end
meta.IOP{:} = '{';
for m=1:length(S.IOP)
    meta.IOP{:} = [meta.IOP{:} num2str(S.IOP(m)) ','];
end
meta.IOP{:}(end) = '}';

%the following changes with every time frame
meta = getDynamicInterfileFields(meta, S, 1);

%remove unused fields
fn = fieldnames(meta);
for z=1:length(fn)
    if isempty(meta.(fn{z})) || (isstruct(meta.(fn{z})) && isempty(fieldnames(meta.(fn{z}))))
        eval(['meta = rmfield(meta,''' fn{z} ''');']);
    end
end

for n=1:nFrames
    
    imgFile = imageFile;
    if nFrames>1
        temp = regexp(imgFile,'\.','split');
        temp{end-1} = [temp{end-1} '_' num2str(n)];
        try
            imgFile = strjoin(temp,'.');
        catch
            for m=1:numel(temp)
                if m==1
                    imgFile = temp{m};
                else
                    imgFile = [imgFile '.' temp{m}];
                end
            end
        end
    end
    
    if exist('h','var') && ishandle(h) && strcmp(get(h,'type'),'figure')
        waitbar(double(n/nFrames),h)
    end
    
    %write image to binary file
    fid = fopen(imgFile,'w');
    %note the flip, this will be wrong if input is not from 'readImages'
    if isfield(S,'volume')
        try
            fwrite(fid, flip(S.volume(:,:,:,n),3), 'single');
        catch
            fwrite(fid, flipdim(S.volume(:,:,:,n),3), 'single');
        end
    elseif isfield(S,'volumes')
        try
            fwrite(fid, flip(S.volumes(:,:,:,n),3), 'single');
        catch
            fwrite(fid, flipdim(S.volumes(:,:,:,n),3), 'single');
        end
    end
    fclose(fid);
    
    %write information to header
    meta = getDynamicInterfileFields(meta, S, n);
    hdrFile = [imgFile '.hdr'];
    writeInterfileHeader(hdrFile, meta);
    
end

if exist('h','var') && ishandle(h) && strcmp(get(h,'type'),'figure')
    close(h)
end
drawnow

end

function [S, imageFile] = parseInputs(varargin)

%check what core architecture for file path character
test = computer;
if strcmpi(test(1:2),'PC')
    pathChar = '\';
else
    pathChar = '/';
end

for n=1:nargin
    if ischar(varargin{n})
        imageFile = varargin{n};
    elseif isstruct(varargin{n})
        S = varargin{n};
    else
        error('Wrong input format')
    end
end
if ~exist('imageFile','var')
    error('Must define image file name')
end
if ~exist('S','var')
    error('Must input structure (see ''readImages'' function')
end
if exist(imageFile,'file')==7
    if ~strcmpi(imageFile(end),pathChar), imageFile = [imageFile pathChar]; end
    imageFile = [imageFile 'Volume.v'];
end

end

function meta = getStaticInterfileFields(meta, S, varargin)

%this function returns a cell string array of header information

if isfield(S,'gantryModel')
    if strcmpi(S.gantryModel(end-2:end),'mMR') || strcmpi(S.gantryModel,'2008')
        meta.gantryModel{:} = '2008';
    elseif strcmpi(S.gantryModel(end-2:end),'mCT') || strcmpi(S.gantryModel,'1104')
        meta.gantryModel{:} = '1104';
    else
        meta.gantryModel{:} = S.gantryModel;
    end
end
if isfield(S,'tracer')
    meta.tracer{:} = S.tracer;
end
if isfield(S,'patientOrientation')
    meta.patientOrientation{:} = S.patientOrientation;
end
if isfield(S,'imSz1')
    meta.imSz1{:} = num2str(S.imSz1);
elseif isfield(S,'volume')
    meta.imSz1{:} = num2str(size(S.volume,1));
elseif isfield(S,'volumes')
    meta.imSz1{:} = num2str(size(S.volumes,1));
end
if isfield(S,'imSz2')
    meta.imSz2{:} = num2str(S.imSz2);
elseif isfield(S,'volume')
    meta.imSz2{:} = num2str(size(S.volume,2));
elseif isfield(S,'volumes')
    meta.imSz2{:} = num2str(size(S.volumes,2));
end
if isfield(S,'imSz3')
    meta.imSz3{:} = num2str(S.imSz3);
elseif isfield(S,'volume')
    meta.imSz3{:} = num2str(size(S.volume,3));
elseif isfield(S,'volumes')
    meta.imSz3{:} = num2str(size(S.volumes,3));
end
if isfield(S,'bedRange1') && isfield(meta,'imSz1')
    meta.voxDim1{:} = num2str(abs(diff(S.bedRange1))/str2double(meta.imSz1{:}),6);
elseif isfield(S,'FOV1') && isfield(meta,'imSz1')
    meta.voxDim1{:} = num2str(abs(S.FOV1)/str2double(meta.imSz1{:}));
end
if isfield(S,'bedRange2') && isfield(meta,'imSz2')
    meta.voxDim2{:} = num2str(abs(diff(S.bedRange2))/str2double(meta.imSz2{:}),6);
elseif isfield(S,'FOV2') && isfield(meta,'imSz2')
    meta.voxDim2{:} = num2str(abs(S.FOV2)/str2double(meta.imSz2{:}));
end
if isfield(S,'bedRange3') && isfield(meta,'imSz3')
    meta.voxDim3{:} = num2str(abs(diff(S.bedRange3))/str2double(meta.imSz3{:}),6);
elseif isfield(S,'FOV3') && isfield(meta,'imSz3')
    meta.voxDim3{:} = num2str(abs(S.FOV3)/str2double(meta.imSz3{:}));
end
if isfield(S,'injectedDose')
    meta.injectedDose{:} = S.injectedDose;
end
if isfield(S,'units')
    meta.units{:} = S.units;
    if strcmpi(meta.units{:},'bqml')
        meta.units{:} = 'Bq/ml';
    elseif strcmpi(meta.units{:},'propcps')
        meta.units{:} = 'cnts/sec';
    elseif strcmpi(meta.units{:},'1/cm')
        meta.units{:} = '1/cm';
    else
        error('Need to define quantification units string.')
    end
end

%translate bedRanges along patient axes (from 'readImages') to scanner axes
[~, bedRange2, bedRange3] = getScannerBedRanges(...
    S.IOP, S.IPP, str2double(meta.imSz1{:}), str2double(meta.imSz2{:}), ...
    str2double(meta.voxDim1{:}), str2double(meta.voxDim2{:}));

%get scanner axes voxel dimensions
[~, ~, voxDimZ] = rotateValues(...
    str2double(meta.voxDim1{:}),str2double(meta.voxDim2{:}),str2double(meta.voxDim3{:}),S.IOP);

if isfield(S,'gantryModel') && ...
        (strcmpi(S.gantryModel(end-2:end),'mCT') || strcmpi(S.gantryModel,'1104'))
    %need to swap to interfile convention
    try
        S.bedRange3 = flip(S.bedRange3);
        bedRange3 = flip(bedRange3);
    catch
        S.bedRange3 = fliplr(S.bedRange3);
        bedRange3 = fliplr(bedRange3);
    end
    %offsets
    [~,y] =  getGantryOffsets(S.gantryModel, S.patientOrientation, S.IOP, str2double(meta.voxDim1{:}), str2double(meta.voxDim1{:}));
    meta.yOffset{:} = num2str(y);
    %bed starts and ends
    if bedRange3(1)<bedRange3(2)
        meta.hBedStart{:} = num2str(bedRange3(1)+abs(voxDimZ)/2);
        meta.hBedEnd{:} = num2str(bedRange3(2)-abs(voxDimZ)/2);
    else
        meta.hBedStart{:} = num2str(bedRange3(1)-abs(voxDimZ)/2);
        meta.hBedEnd{:} = num2str(bedRange3(2)+abs(voxDimZ)/2);
    end
    meta.vBedStart{:} = num2str(-mean(bedRange2));
    %account for offset in actual bed boundaries because it's not read from header
    if isfield(meta,'yOffset')
        meta.vBedStart{:} = num2str(str2double(meta.vBedStart{:})+str2double(meta.yOffset{:}));
    end
    if isfield(meta,'zOffset')
        meta.hBedStart{:} = num2str(str2double(meta.hBedStart{:})+str2double(meta.zOffset{:}));
        meta.hBedEnd{:} = num2str(str2double(meta.hBedEnd{:})+str2double(meta.zOffset{:}));
    end
elseif isfield(S,'gantryModel') && ...
        (strcmpi(S.gantryModel(end-2:end),'mMR') || strcmpi(S.gantryModel,'2008'))
    %offsets
    [x,y,z] =  getGantryOffsets(S.gantryModel, S.patientOrientation, S.IOP, str2double(meta.voxDim1{:}), str2double(meta.voxDim1{:}));
    meta.xOffset{:} = num2str(x);
    meta.yOffset{:} = num2str(y);
    meta.zOffset{:} = num2str(z);
    bedRange2 = bedRange2-y;
    S.bedRange2 = S.bedRange2-y;
    if S.locZ(1)<S.locZ(end) && S.bedRange3(1)>S.bedRange3(end) || ...
            S.locZ(1)>S.locZ(end) && S.bedRange3(1)<S.bedRange3(end)
        bedRange3 = bedRange3+z;
    else
        bedRange3 = bedRange3-z;
    end
    S.bedRange3 = S.bedRange3-z;
    %mMR requires a calculation to translate actual bed range to header fields
    temp = cross(S.IOP(1:3),S.IOP(4:6));
    S.bedRange3 = (temp(3) * S.bedRange3) - 1 * ((127 - 1) / 2.0 * 2.03125);
    bedRange3 = (temp(3) * bedRange3) - 1 * ((127 - 1) / 2.0 * 2.03125);
    %this accounts for volume flip
    if temp(3)*bedRange3(end)>temp(3)*bedRange3(1)
        try
            S.bedRange3 = flip(S.bedRange3);
            bedRange3 = flip(bedRange3);
        catch
            S.bedRange3 = fliplr(S.bedRange3);
            bedRange3 = fliplr(bedRange3);
        end
    end
    
    %bed starts and ends
    if bedRange3(1)<bedRange3(2)
        meta.hBedStart{:} = num2str(bedRange3(1)+abs(voxDimZ)/2);
        meta.hBedEnd{:} = num2str(bedRange3(2)-abs(voxDimZ)/2);
    else
        meta.hBedStart{:} = num2str(bedRange3(1)-abs(voxDimZ)/2);
        meta.hBedEnd{:} = num2str(bedRange3(2)+abs(voxDimZ)/2);
    end
    meta.vBedStart{:} = num2str(-mean(bedRange2));
else %need to swap to interfile convention (this was done for mMR by * -1)
    try
        S.bedRange3 = flip(S.bedRange3);
        bedRange3 = flip(bedRange3);
    catch
        S.bedRange3 = fliplr(S.bedRange3);
        bedRange3 = fliplr(bedRange3);
    end
    %bed starts and ends
    if bedRange3(1)<bedRange3(2)
        meta.hBedStart{:} = num2str(bedRange3(1)+abs(voxDimZ)/2);
        meta.hBedEnd{:} = num2str(bedRange3(2)-abs(voxDimZ)/2);
    else
        meta.hBedStart{:} = num2str(bedRange3(1)-abs(voxDimZ)/2);
        meta.hBedEnd{:} = num2str(bedRange3(2)+abs(voxDimZ)/2);
    end
    meta.vBedStart{:} = num2str(-mean(bedRange2));
    %account for offset in actual bed boundaries because it's not read from header
    if isfield(meta,'yOffset')
        meta.vBedStart{:} = num2str(str2double(meta.vBedStart{:})+str2double(meta.yOffset{:}));
    end
    if isfield(meta,'zOffset')
        meta.hBedStart{:} = num2str(str2double(meta.hBedStart{:})+str2double(meta.zOffset{:}));
        meta.hBedEnd{:} = num2str(str2double(meta.hBedEnd{:})+str2double(meta.zOffset{:}));
    end
end

if isfield(S,'calibrationFctr')
    meta.calibrationFctr{:} = num2str(S.calibrationFctr,6);
end
if isfield(S,'studyDate')
    meta.studyDate{:} = S.studyDate;
    test = regexp(meta.studyDate{:},':','split');
    if length(test)==1
        meta.studyDate{:} = [meta.studyDate{:}(1:4) ':' meta.studyDate{:}(5:6) ':' meta.studyDate{:}(7:8)];
    end
end
if isfield(S,'studyTime')
    meta.studyTime{:} = S.studyTime;
    test = regexp(meta.studyTime{:},':','split');
    if length(test)==1
        meta.studyTime{:} = [meta.studyTime{:}(1:2) ':' meta.studyTime{:}(3:4) ':' meta.studyTime{:}(5:6)];
    end
end
if isfield(S,'injectionDate')
    meta.injDate{:} = S.injectionDate;
    test = regexp(meta.injDate{:},':','split');
    if length(test)==1
        meta.injDate{:} = [meta.injDate{:}(1:4) ':' meta.injDate{:}(5:6) ':' meta.injDate{:}(7:8)];
    end
end
if isfield(S,'injectionTime')
    meta.injTime{:} = S.injectionTime;
    test = regexp(meta.injTime{:},':','split');
    if length(test)==1
        meta.injTime{:} = [meta.injTime{:}(1:2) ':' meta.injTime{:}(3:4) ':' meta.injTime{:}(5:6)];
    end
end

if ~isfield(meta,'gantryModel'), meta.gantryModel{:} = ''; end
if ~isfield(meta,'tracer'), meta.tracer{:} = ''; end
if ~isfield(meta,'patientOrientation'), meta.patientOrientation{:} = ''; end
if ~isfield(meta,'IOP'), meta.IOP{:} = ''; end
if ~isfield(meta,'imSz1'), meta.imSz1{:} = ''; end
if ~isfield(meta,'imSz2'), meta.imSz2{:} = ''; end
if ~isfield(meta,'imSz3'), meta.imSz3{:} = ''; end
if ~isfield(meta,'voxDim1'), meta.voxDim1{:} = ''; end
if ~isfield(meta,'voxDim2'), meta.voxDim2{:} = ''; end
if ~isfield(meta,'voxDim3'), meta.voxDim3{:} = ''; end
if ~isfield(meta,'injectedDose'), meta.injectedDose{:} = ''; end
if ~isfield(meta,'units'), meta.units{:} = ''; end
if ~isfield(meta,'hBedStart'), meta.hBedStart{:} = ''; end
if ~isfield(meta,'hBedEnd'), meta.hBedEnd{:} = ''; end
if ~isfield(meta,'vBedStart'), meta.vBedStart{:} = ''; end
if ~isfield(meta,'calibrationFctr'), meta.calibrationFctr{:} = ''; end
if ~isfield(meta,'startTime'), meta.startTime{:} = ''; end
if ~isfield(meta,'endTime'), meta.endTime{:} = ''; end
if ~isfield(meta,'prompts'), meta.prompts{:} = ''; end
if ~isfield(meta,'randoms'), meta.randoms{:} = ''; end
if ~isfield(meta,'studyDate'), meta.studyDate{:} = ''; end
if ~isfield(meta,'studyTime'), meta.studyTime{:} = ''; end
if ~isfield(meta,'injDate'), meta.injDate{:} = ''; end
if ~isfield(meta,'injTime'), meta.injTime{:} = ''; end

end

function meta = getDynamicInterfileFields(meta, S, n)

if isfield(S,'startTime')
    meta.startTime{:} = num2str(S.startTime(n));
elseif isfield(S,'startTimes')
    meta.startTime{:} = num2str(S.startTimes(n));
end
if isfield(S,'endTime')
    meta.endTime{:} = num2str(S.endTime(n));
elseif isfield(S,'endTimes')
    meta.endTime{:} = num2str(S.endTimes(n));
end
if isfield(S,'prompts')
    meta.prompts{:} = num2str(S.prompts(n));
end
if isfield(S,'randoms')
    meta.randoms{:} = num2str(S.randoms(n));
end
if isfield(S,'volume')
    meta.maxPixel{:} = sprintf('%.1f',max(max(max(S.volume(:,:,:,n)))));
    meta.minPixel{:} = sprintf('%.1f',min(min(min(S.volume(:,:,:,n)))));
elseif isfield(S,'volume')
    meta.maxPixel{:} = sprintf('%.1f',max(max(max(S.volumes(:,:,:,n)))));
    meta.minPixel{:} = sprintf('%.1f',min(min(min(S.volumes(:,:,:,n)))));
end

end

function out = initiateInterfileStruct

out = struct;
out.gantryModel = {''};
out.tracer = {''};
out.patientOrientation = {''};
out.IOP = {''};
out.imSz1 = {''};
out.imSz2 = {''};
out.imSz3 = {''};
out.voxDim1 = {''};
out.voxDim2 = {''};
out.voxDim3 = {''};
out.injectedDose = {''};
out.units = {''};
out.hBedStart = {''};
out.hBedEnd = {''};
out.vBedStart = {''};
out.calibrationFctr = {''};
out.startTime = {''};
out.endTime = {''};
out.prompts = {''};
out.randoms = {''};
out.studyDate = {''};
out.studyTime = {''};
out.injDate = {''};
out.injTime = {''};
out.maxPixel = {''};
out.minPixel = {''};

end

function writeInterfileHeader(hdrFile, meta)

%inputs are header file name and structure of cell string values

%check what core architecture for file path character
test = computer;
if strcmpi(test(1:2),'PC')
    pathChar = '\';
else
    pathChar = '/';
end

if ~strcmpi(hdrFile(end-3:end),'.hdr'), hdrFile = [hdrFile '.hdr']; end

temp = regexp(hdrFile,pathChar,'split');
temp = temp{end};
temp = regexp(temp,'\.','split');
try
    imgFile{:} = strjoin(temp(1:end-1),'.');
catch
    for n=1:numel(temp)-1
        if n==1
            imgFile{:} = temp{n};
        else
            imgFile{:} = [imgFile{:} '.' temp{n}];
        end
    end
end

if ~isempty(meta.prompts{:}) && ~isempty(meta.randoms{:})
    trues = num2str(str2double(meta.prompts{:})-str2double(meta.randoms{:}));
else
    trues = '';
end

fid = fopen(hdrFile,'w');

%write header file information
fprintf(fid,'!INTERFILE:=\r\n');
fprintf(fid,'%%comment:=\r\n');
fprintf(fid,['!originating system:=' meta.gantryModel{:} '\r\n']);
fprintf(fid,'%%SMS-MI header name space:=\r\n');
fprintf(fid,'%%SMS-MI version number:=\r\n');
fprintf(fid,'\r\n');
fprintf(fid,'!GENERAL DATA:=\r\n');
fprintf(fid,'%%sinogram header file:=\r\n');
fprintf(fid,'%%sinogram data file:=\r\n');
fprintf(fid,['!name of data file:=' imgFile{:} '\r\n']);
fprintf(fid,'\r\n');
fprintf(fid,'!GENERAL IMAGE DATA:=\r\n');
fprintf(fid,['%%study date (yyyy:mm:dd):=' meta.studyDate{:} '\r\n']);
fprintf(fid,['%%study time (hh:mm:ss GMT+00:00):=' meta.studyTime{:} '\r\n']);
fprintf(fid,'isotope name:=\r\n');
fprintf(fid,'isotope gamma halflife (sec):=\r\n');
fprintf(fid,'isotope branching factor:=\r\n');
fprintf(fid,['radiopharmaceutical:=' meta.tracer{:} '\r\n']);
fprintf(fid,['%%tracer injection date (yyyy:mm:dd):=' meta.injDate{:} '\r\n']);
fprintf(fid,['%%tracer injection time (hh:mm:ss GMT+00:00):=' meta.injTime{:} '\r\n']);
fprintf(fid,['tracer activity at time of injection (Bq):=' num2str(meta.injectedDose{:},5) '\r\n']);
fprintf(fid,'relative time of tracer injection (sec):=\r\n');
fprintf(fid,'injected volume (ml):=\r\n');
fprintf(fid,'image data byte order:=LITTLEENDIAN\r\n');
fprintf(fid,['%%patient orientation:=' meta.patientOrientation{:} '\r\n']);
fprintf(fid,['%%image orientation:=' meta.IOP{:} '\r\n']);
fprintf(fid,'!PET data type:=image\r\n');
fprintf(fid,'number format:=float\r\n');
fprintf(fid,'!number of bytes per pixel:=4\r\n');
fprintf(fid,'number of dimensions:=3\r\n');
fprintf(fid,'matrix axis label[1]:=x\r\n');
fprintf(fid,'matrix axis label[2]:=y\r\n');
fprintf(fid,'matrix axis label[3]:=z\r\n');
fprintf(fid,['matrix size[1]:=' meta.imSz1{:} '\r\n']);
fprintf(fid,['matrix size[2]:=' meta.imSz2{:} '\r\n']);
fprintf(fid,['matrix size[3]:=' meta.imSz3{:} '\r\n']);
fprintf(fid,['scale factor (mm/pixel) [1]:=' meta.voxDim1{:} '\r\n']);
fprintf(fid,['scale factor (mm/pixel) [2]:=' meta.voxDim2{:} '\r\n']);
fprintf(fid,['scale factor (mm/pixel) [3]:=' meta.voxDim3{:} '\r\n']);
fprintf(fid,'horizontal bed translation:=stepped\r\n');
fprintf(fid,['start horizontal bed position (mm):=' sprintf('%.3f',str2double(meta.hBedStart{:})) '\r\n']);
fprintf(fid,['end horizontal bed position (mm):=' sprintf('%.3f',str2double(meta.hBedEnd{:})) '\r\n']);
fprintf(fid,['start vertical bed position (mm):=' sprintf('%.3f',str2double(meta.vBedStart{:})) '\r\n']);
fprintf(fid,'%%reconstruction diameter (mm):=\r\n');
fprintf(fid,['quantification units:=' meta.units{:} '\r\n']);
fprintf(fid,['%%scanner quantification factor (Bq*s/ECAT counts):=' meta.calibrationFctr{:} '\r\n']);
fprintf(fid,'%%decay correction:=\r\n');
fprintf(fid,['%%decay correction reference date (yyyy:mm:dd):=' meta.studyDate{:} '\r\n']);
fprintf(fid,['%%decay correction reference time (hh:mm:ss GMT+00:00):=' meta.studyTime{:} '\r\n']);
fprintf(fid,'slice orientation:=\r\n');
fprintf(fid,'method of reconstruction:=\r\n');
fprintf(fid,'%%PSF axial sigma (mm):=\r\n');
fprintf(fid,'%%PSF axial cutoff:=\r\n');
fprintf(fid,'%%number of radial PSF coefficients:=\r\n');
fprintf(fid,'%%PSF radial left coefficient (bins) [1]:=\r\n');
fprintf(fid,'%%PSF radial left coefficient (1/mm) [2]:=\r\n');
fprintf(fid,'%%PSF radial left coefficient (1/mm^2) [3]:=\r\n');
fprintf(fid,'%%PSF radial left coefficient (1/mm^3) [4]:=\r\n');
fprintf(fid,'%%PSF radial left coefficient (1/mm^4) [5]:=\r\n');
fprintf(fid,'%%PSF radial right coefficient (bins) [1]:=\r\n');
fprintf(fid,'%%PSF radial right coefficient (1/mm) [2]:=\r\n');
fprintf(fid,'%%PSF radial right coefficient (1/mm^2) [3]:=\r\n');
fprintf(fid,'%%PSF radial right coefficient (1/mm^3) [4]:=\r\n');
fprintf(fid,'%%PSF radial right coefficient (1/mm^4) [5]:=\r\n');
fprintf(fid,'%%PSF radial cutoff:=\r\n');
if isfield(meta,'xOffset')
    fprintf(fid,['%%gantry offset (mm) [1]:=' meta.xOffset{:} '\r\n']);
else
    fprintf(fid,'%%gantry offset (mm) [1]:=0\r\n');
end
if isfield(meta,'yOffset')
    fprintf(fid,['%%gantry offset (mm) [2]:=' meta.yOffset{:} '\r\n']);
else
    fprintf(fid,'%%gantry offset (mm) [2]:=0\r\n');
end
if isfield(meta,'zOffset')
    fprintf(fid,['%%gantry offset (mm) [3]:=' meta.zOffset{:} '\r\n']);
else
    fprintf(fid,'%%gantry offset (mm) [3]:=0\r\n');
end
fprintf(fid,'%%gantry offset pitch (degrees):=\r\n');
fprintf(fid,'%%gantry offset yaw (degrees):=\r\n');
fprintf(fid,'%%gantry offset roll (degrees):=\r\n');
fprintf(fid,'number of iterations:=\r\n');
fprintf(fid,'%%number of subsets:=\r\n');
fprintf(fid,'filter name:=\r\n');
fprintf(fid,'%%xy-filter (mm):=\r\n');
fprintf(fid,'%%z-filter (mm):=\r\n');
fprintf(fid,'%%filter order:=\r\n');
fprintf(fid,'%%image zoom:=\r\n');
fprintf(fid,'%%x-offset (mm):=\r\n');
fprintf(fid,'%%y-offset (mm):=\r\n');
fprintf(fid,'applied corrections:=\r\n');
fprintf(fid,'method of attenuation correction:=\r\n');
fprintf(fid,'%%CT coverage:=\r\n');
fprintf(fid,'method of scatter correction:=\r\n');
fprintf(fid,'%%method of random correction:=\r\n');
fprintf(fid,'%%TOF mashing factor:=\r\n');
fprintf(fid,'number of energy windows:=\r\n');
fprintf(fid,'%%energy window lower level (keV) [1]:=\r\n');
fprintf(fid,'%%energy window upper level (keV) [1]:=\r\n');
fprintf(fid,'%%coincidence window width (ns):=\r\n');
fprintf(fid,'\r\n');
fprintf(fid,'!IMAGE DATA DESCRIPTION:=\r\n');
fprintf(fid,'!total number of data sets:=1\r\n');
fprintf(fid,['!image duration (sec):=' ...
    num2str(str2double(meta.endTime{:})-str2double(meta.startTime{:})) '\r\n']);
fprintf(fid,['!image relative start time (sec):=' meta.startTime{:} '\r\n']);
fprintf(fid,['total prompts:=' meta.prompts{:} '\r\n']);
fprintf(fid,['%%total randoms:=' meta.randoms{:} '\r\n']);
fprintf(fid,['%%total net trues:=' trues '\r\n']);
fprintf(fid,'%%GIM loss fraction:=\r\n');
fprintf(fid,'%%PDR loss fraction:=\r\n');
fprintf(fid,'%%total uncorrected singles rate:=\r\n');
fprintf(fid,'%%image slope:=\r\n');
fprintf(fid,'%%image intercept:=\r\n');
fprintf(fid,['maximum pixel count:=' meta.maxPixel{:} '\r\n']);
fprintf(fid,['minimum pixel count:=' meta.minPixel{:} '\r\n']);
fprintf(fid,'\r\n');
fprintf(fid,'%%SUPPLEMENTARY ATTRIBUTES:=\r\n');

fclose(fid);

end