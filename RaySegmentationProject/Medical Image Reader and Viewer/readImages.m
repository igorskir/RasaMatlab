function out = readImages(varargin)

% this function reads medical imaging data from an image file or folder
% important information is stored in a convenient structure for MATLAB
% this structure is used as the input for many other functions
%
% function can be called with no input argumets, will prompt for file selection
% or a specific image file can be passed
% dynamic series can be loaded by passing the series directory path (no file)
%
% output volume X, Y, and Z axes count columns, rows and slices, respectively
% this is intuitive, but is rotated 90 degrees relative to DICOM convention
%
% NOTE: output ranges are in image coordinates
%      which are not scanner coordinates (unless IOP==[1; 0; 0; 0; 1; 0])

warning('off','MATLAB:warn_r14_stucture_assignment')

out = struct;

doAlign = 0;

inStruct = parseInputs(varargin{:});
if isfield(inStruct,'align')
    doAlign = 1;
end

if isfield(inStruct,'imageFiles') && isempty(inStruct.imageFiles)
    out = [];  %default return if no image(s) loaded
    return;
end

imageDir = inStruct.imageDir;
if isfield(inStruct,'imageFiles'), imageFiles = inStruct.imageFiles; end
fileType = inStruct.fileType;

runFlag = 1;
while runFlag
    try
        if ~isempty(strfind(fileType,'dcm')) || ~isempty(strfind(fileType,'ima'))
            if exist('imageFiles','var')
                dicomOut = processDicom(imageDir, imageFiles);
            else
                dicomOut = processDicom(imageDir);
            end
        end
        
        runFlag = 0;
        
    catch
        %         disp(' ')
        close(findall(0,'Tag','TMWWaitbar'))
        disp('File read error, try again with another file')
        %clean up
        clear volumes
        if exist('h','var') && ishandle(h) && strcmp(get(h,'type'),'figure')
            close(h)
        end
        %try a new file
        filelist = {'*.img;*.v;*.dcm;*.ima', ...
            'All Image Files (*.img, *.v, *.dcm, *.ima)'; ...
            '*.*','All Files (*.*)'};
        [imageFile, imageDir] = uigetfile(filelist,'Select Image File',imageDir);
        if imageFile==0
            out = [];  %default return if no image(s) loaded
            return;
        end
        varargin = {[imageDir imageFile]};
        
        inStruct = parseInputs(varargin{:});
        if isempty(inStruct.imageFiles)
            out = [];  %default return if no image(s) loaded
            return;
        end
        
        imageDir = inStruct.imageDir;
        imageFiles = inStruct.imageFiles;
        fileType = inStruct.fileType;
        
        runFlag = 1;
    end
end

%% this starts the code for all images

nSeries = 0;
%find all returned image structures
if exist('dicomOut','var')
    nSeries = nSeries+dicomOut.nSeries;
end

out.nSeries = nSeries;

if exist('dicomOut','var')
    fn = fieldnames(dicomOut);
    fn(strcmpi(fn,'nSeries')) = [];
    for n=1:length(fn)
        %check if series name already exists in out structure
        temp = fn{n};
        if ismember(temp,fieldnames(out))
            temp = [fn{n} '_dicom'];
        end
        out.(temp) = dicomOut.(fn{n});
    end
end
sNames = fieldnames(out);
sNames(strcmpi(sNames,'nSeries')) = [];

% h = waitbar(0,'Mapping voxel locations...');

for n=1:nSeries
    S = out.(sNames{n});    
    if isfield(S,'volumes') %only needed for images, i.e. not spectroscopy        
        %assure consistently spaced slices
        vect3 = cross(S.IOP(1:3),S.IOP(4:6)); %image slice axis
        temp = vect3*S.voxDim3; %find image slice spacing contribution to real scanner axes
        temp1 = temp(1); %voxDim3 contribution along axis 1
        temp2 = temp(2); %voxDim3 contribution along axis 2
        temp3 = temp(3); %voxDim3 contribution along axis 3
        if temp'*S.IPP(:,1)<temp'*S.IPP(:,end)
            if size(S.IPP,2)>1
                idx = round((S.IPP-repmat(S.IPP(:,1),1,S.imSz3))./repmat([temp1; temp2; temp3],1,S.imSz3));
                realIdx = find(~isnan(idx(:,round((size(idx,2)+1)/2))) & ~isinf(idx(:,round((size(idx,2)+1)/2))) & idx(:,round((size(idx,2)+1)/2))~=0);
                idx = idx(realIdx(1),:);
            else
                idx = 0;
            end
            S.IPP = repmat(S.IPP(:,1),1,S.imSz3)+repmat(idx,3,1).*repmat([temp1; temp2; temp3],1,S.imSz3);
        else
            if size(S.IPP,2)>1
                idx = round((S.IPP-repmat(S.IPP(:,end),1,S.imSz3))./repmat([temp1; temp2; temp3],1,S.imSz3));
                realIdx = find(~isnan(idx(:,round(size(idx,2)/2))) & ~isinf(idx(:,round(size(idx,2)/2))) & idx(:,round(size(idx,2)/2))~=0);
                idx = idx(realIdx(1),:);
            else
                idx = 0;
            end
            S.IPP = repmat(S.IPP(:,end),1,S.imSz3)+repmat(idx,3,1).*repmat([temp1; temp2; temp3],1,S.imSz3);
        end
        
        %hack, ECAT series would break here because no IPP or IOP
        if ~isfield(S,'IOP'), S.IOP = [1; 0; 0; 0; 1; 0]; end
        %bed ranges may NOT be correct until this step
        [S.bedRange1, S.bedRange2, S.bedRange3, S.locX, S.locY, S.locZ] = ...
            getVolumeBedRanges(S.IOP, S.IPP, S.imSz1, S.imSz2, S.voxDim1, S.voxDim2, S.voxDim3);
        if doAlign
            if size(S.volumes,3)>1
                if ~all(round(S.IOP*1e4)/1e4==[1;0;0;0;1;0]) || ...
                        S.bedRange1(1)>S.bedRange1(2) || ...
                        S.bedRange2(1)>S.bedRange2(2) || ...
                        S.bedRange3(1)>S.bedRange3(2)
                    disp('Aligning volume...')
                    %create structure for 'alignVolume' function
                    tempStruct = struct;
                    tempStruct.volumes = S.volumes;
                    tempStruct.IOP = S.IOP;
                    tempStruct.IPP = S.IPP;
                    tempStruct.bedRange1 = S.bedRange1;
                    tempStruct.bedRange2 = S.bedRange2;
                    tempStruct.bedRange3 = S.bedRange3;
                    %do the alignment
                    temp = alignVolume(tempStruct);
                    %update variables
                    S.volumes = temp.volumes;
                    S.IOP = temp.IOP;
                    S.IPP = temp.IPP;
                    S.locX = temp.locX;
                    S.locY = temp.locY;
                    S.locZ = temp.locZ;
                    S.bedRange1 = temp.bedRange1;
                    S.bedRange2 = temp.bedRange2;
                    S.bedRange3 = temp.bedRange3;
                end
            else
                disp('Slices will not be aligned')
            end
        end
        %convert voxel locations to single to conserve memory
        if isfield(S,'locX') && ~isa(S.locX,'single')
            S.locX = single(S.locX);
            S.locY = single(S.locY);
            S.locZ = single(S.locZ);
        end
        out.(sNames{n}) = rmfield(S,{'voxDim1','voxDim2','voxDim3'});
    end    
    %     if exist('h','var') && ishandle(h) && strcmp(get(h,'type'),'figure')
    %         waitbar(double(n/nSeries),h)
    %     end
end

% if exist('h','var') && ishandle(h) && strcmp(get(h,'type'),'figure')
%     close(h)
%     drawnow
% end

%now check if only 1 total complete series found, if so move it to output structure "root"
fn = fieldnames(out);
fn(strcmpi(fn,'nSeries')) = [];
if out.nSeries==1 && length(fn)==1
    sName = fn{1};
    fn = fieldnames(out.(sName));
    for m=1:length(fn)
        out.(fn{m}) = out.(sName).(fn{m});
    end
    %remove series sub structure
    eval(['out = rmfield(out,''' sName ''');'])
end

end

function out = parseInputs(varargin)

%check what core architecture for file path character
test = computer;
if strcmpi(test(1:2),'PC')
    pathChar = '\';
else
    pathChar = '/';
end

if any(strcmpi(varargin,'align'))
    out.align = 1;
    varargin(strcmpi(varargin,'align')) = [];
end
in = varargin;

useDir = 0;
if ~isempty(varargin) && exist(varargin{1},'file')==7
    useDir = 1;
end

out.imageFiles = [];

filelist = {'*.img;*.v;*.dcm;*.ima', ...
    'All Image Files (*.img, *.v, *.dcm, *.ima)'; ...
    '*.*','All Files (*.*)'};

%if no input, prompt user to select image directory
if ~exist('in','var') || isempty(in)
    [imageFile, imageDir] = uigetfile(filelist,'Select Image File',cd);
    if imageFile==0, return; end
    in = [imageDir imageFile];
    useDir = 0;
else
    if iscell(in), in = in{:}; end
end

if ~(exist(in, 'file')==2) && ~(exist(in, 'file')==7)
    error('Input argument is not a valid file or path.')
end

while ~exist('imageFiles','var')
    if exist(in, 'file')==2
        %input is a specific image file
        imageDir = regexp(in,pathChar,'split');
        imageFile = imageDir{end};
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
        if strcmpi(imageDir,'')
            imageDir = [cd pathChar];
        elseif ~strcmpi(imageDir(end),pathChar)
            imageDir = [imageDir pathChar];
        end
        imageFiles = {imageFile};
        useDir = 0;
    else
        imageDir = in;
        if ~strcmpi(imageDir(end),pathChar), imageDir = [imageDir pathChar]; end
        allFiles = dir([imageDir '*.*']);
        allFiles([allFiles.isdir]) = [];
        allFiles = {allFiles.name};
        %find extensions
        temp = regexp(allFiles,'\.','split');
        ext = cell(length(temp),1);
        for n=1:length(temp)
            ext(n) = temp{n}(end);
        end
        
        ecatIdx = strcmpi(ext,'img');
        interfileIdx = strcmpi(ext,'v');
        dicomIdx = strcmpi(ext,'dcm') | strcmpi(ext,'ima');
        otherIdx = ~(ecatIdx | interfileIdx | dicomIdx);
        for n=1:sum(otherIdx)
            idx = find(otherIdx);
            test = dicm_hdr([imageDir allFiles{idx(n)}]);
            if ~isempty(test), dicomIdx(idx(n)) = 1; end
        end
        %construct string that holds all extension suffixes
        extString = '';
        if sum(ecatIdx)>0
            extString = [extString 'img'];
        end
        if sum(interfileIdx)>0
            extString = [extString 'v'];
        end
        if sum(dicomIdx)>0
            extString = [extString 'dcm'];
        end
        if isempty(extString)
            %no image files found in input directory
            [imageFile, imageDir] = uigetfile(filelist,'Select Image File',imageDir);
            if imageFile==0, return; end
            imageFiles = {imageFile};
            useDir = 0;
        else
            %dummy variable to satisfy 'while'
            imageFiles = [];
            useDir = 1;
        end
    end
    try
        imageFiles = validateImageFiles(imageDir, imageFiles);
    catch
        %prompt for new image file
        [imageFile, imageDir] = uigetfile(filelist,'Select Image File',imageDir);
        if imageFile==0, return; end
        in = [imageDir imageFile];
    end
end

if useDir == 0;
    %sort files
    temp = regexp(imageFiles,'[0123456789]','match');
    num = zeros(length(imageFiles),1);
    for n=1:length(imageFiles)
        num(n) = str2double([temp{n}{:}]);
    end
    [~,idx] = sort(num);
    imageFiles = imageFiles(idx);
    %get filetype
    temp = [imageDir imageFiles{1}];
    temp = regexp(temp,'\.','split');
    ext = lower(temp{end});
    out.imageFiles = imageFiles;
    out.fileType = ext;
    %in case extensionless dicom
    test = dicm_hdr([imageDir imageFiles{1}]);
    if ~isempty(test), out.fileType = 'dcm'; end
else
    out.fileType = extString;
end
out.imageDir = imageDir;

if useDir %if passed a directory first
    out = rmfield(out,'imageFiles');
end

end

function imageFiles = validateImageFiles(imageDir, imageFiles)

for n=1:length(imageFiles)
    
    imageFile = imageFiles{n};
    
    imFileParts = regexp(imageFile,'\.','split');
    ext = imFileParts{end};
    
    if ~strcmpi(ext,'dcm') && ~strcmpi(ext,'ima') && ~strcmpi(ext,'img') && ~strcmpi(ext,'v') && ~strcmpi(ext,'hdr')
        %first see if it is extensionless DICOM file
        test = dicm_hdr([imageDir imageFile]);
        if isempty(test)
            disp(['Unknown image file type: ''' imageFile ''''])
            error(['Unknown image file type: ''' imageFile ''''])
        end
    end
    if strcmpi(ext,'img') || strcmpi(ext,'v') || strcmpi(ext,'hdr')
        %check if input file is header
        if strcmpi(ext,'hdr')
            hdrFile = imageFile;
            try
                imageFile = strjoin(imFileParts(1:end-1),'.');
            catch
                temp = imageFile;
                for m=1:numel(temp)-1
                    if m==1
                        imageFile = temp{m};
                    else
                        imageFile = [imageFile '.' temp{m}];
                    end
                end
            end            
            
            %check if image exists
            if ~(exist([imageDir imageFile],'file')==2)
                error(['Image file does not exist for ''' imageDir hdrFile ''''])
            end
        else
            hdrFile = [imageFile '.hdr'];
            %check if header exists
            if ~(exist([imageDir hdrFile],'file')==2)
                error(['Header file does not exist for ''' imageDir imageFile ''''])
            end
        end
    end
    
    imageFiles{n} = imageFile;
    
end

end

function out = processDicom(imageDir, varargin)

if ~isempty(varargin), imageFiles = varargin{:}; end

if exist('imageFiles','var')
    if iscell(imageFiles), imageFiles = imageFiles{:}; end
    DICOMstruct = readDicomImages([imageDir imageFiles]);
else
    DICOMstruct = readDicomImages(imageDir);
end

nSeries = DICOMstruct.nSeries;
DICOMstruct = rmfield(DICOMstruct,'nSeries');
out.nSeries = nSeries;
sNames = fieldnames(DICOMstruct);
for n=1:length(sNames)
    S = DICOMstruct.(sNames{n});
    if isfield(S,'volume') || isfield(S,'volumes')
        if isfield(S,'volume')
            S.volumes = permute(S.volume,[2 1 3 4]);
            S = rmfield(S,'volume');
        else
            S.volumes = permute(S.volumes,[2 1 3 4]);
        end
        %account for flip in matrix and voxel dimensions
        S.imSz1 = size(S.volumes,1);
        S.imSz2 = size(S.volumes,2);
        S.imSz3 = size(S.volumes,3);
    end
    S.filePath = imageDir;
    if exist('imageFiles','var')
        S.fileName = imageFiles;
    end
    
    out.(sNames{n}) = S;
end

end