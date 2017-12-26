function writeDicomImage(varargin)

%inputs are structure from 'readImages' function and folder for dicom slices

warning('off','all')

%check core architecture for file path character
test = computer;
if strcmpi(test(1:2),'PC')
    pathChar = '\';
else
    pathChar = '/';
end

[S, imageDir] = parseInputs(varargin{:});
if ~strcmpi(imageDir(end),pathChar), imageDir = [imageDir pathChar]; end
mkdir(imageDir)

%check for previous files
prevFiles = dir(imageDir);
prevFiles(strcmpi({prevFiles.name},'.')) = [];
prevFiles(strcmpi({prevFiles.name},'..')) = [];
%only files
prevFiles([prevFiles.isdir]) = [];
if ~isempty(prevFiles)
    warning('Previous files found in image directory. Moving them to folder ''BAK''')
    mkdir([imageDir pathChar 'BAK'])
    for n=1:length(prevFiles)
        movefile([imageDir pathChar prevFiles(n).name],[imageDir pathChar 'BAK' pathChar prevFiles(n).name])
    end
end

%the rotation and flip orients the images so X counts rows and Y counts columns
if isfield(S,'volume')
%     vol = flip(rot90(S.volume,1));
    vol = permute(S.volume,[2 1 3 4]);
elseif isfield(S,'volumes')
%     vol = flip(rot90(S.volumes,1));
    vol = permute(S.volumes,[2 1 3 4]);
else
    error('Input structure should contain a single volume set.')
end
nFrames = size(vol,4);

%get data modality
SOPClassUID = getSOPClassUID(S);

%initiate dicom structure
meta = initiateDicomStruct;
%set static dicom fields
meta = getStaticDicomFields(meta, S);
%fill in missing info
meta.SOPClassUID = SOPClassUID;
if nFrames>1, meta.NumberOfTimeSlices = nFrames; end
%the following will be the same for the entire series
SOPInstanceUID = dicomuid;
meta.SOPInstanceUID = SOPInstanceUID;
SeriesInstanceUID = dicomuid;
meta.SeriesInstanceUID = SeriesInstanceUID;
%the following will change with every slice
meta.ContentTime = char(datetime('now','Format','HHmmss.SSSSSS'));
StudyInstanceUID = dicomuid;
meta.StudyInstanceUID = StudyInstanceUID;
meta.RescaleIntercept = 0;
meta.RescaleSlope = 1;

%remove unused fields
fn = fieldnames(meta);
for z=1:length(fn)
    if isempty(meta.(fn{z})) || (isstruct(meta.(fn{z})) && isempty(fieldnames(meta.(fn{z}))))
        eval(['meta = rmfield(meta,''' fn{z} ''');']);
    end
end

totalFiles = nFrames * meta.NumberOfSlices;
sliceCounter = 1;

disp(['Writing ' num2str(totalFiles) ' dicom files...'])
h = waitbar(0,['Writing ' num2str(totalFiles) ' dicom files...']);

for n=1:nFrames
    %set dynamic dicom fields
    meta = getDynamicDicomFields(meta, S, n);
    %write slices
    for m=1:meta.NumberOfSlices
        sliceName = ['SLICE.' sprintf('%04d',sliceCounter) '.IMA'];
        
        meta.InstanceNumber = sliceCounter;
        
        meta.ImagePositionPatient = S.IPP(:,m);
        if max(abs(S.IOP))==1 && min(abs(S.IOP))==0
            meta.SliceLocation = S.IPP(3,m);
        end
        
        if exist('h','var') && ishandle(h) && strcmp(get(h,'type'),'figure')
            waitbar(double(sliceCounter/totalFiles),h)
        end
        
        %same for entire series
        meta.SOPClassUID = SOPClassUID;
        meta.SOPInstanceUID = SOPInstanceUID;
        meta.SeriesInstanceUID = SeriesInstanceUID;
        %change with every slice
        meta.ContentTime = char(datetime('now','Format','HHmmss.SSSSSS'));
        %increment instance UID
        temp = regexp(StudyInstanceUID,'\.','split');
        stubLngth = 10;
        offset = 0;
        stub = num2str(str2double(temp{end}(end-(stubLngth-1):end))+1);
        while length(stub)>stubLngth+offset %should almost never happen
            offset = offset+1;
            stub = num2str(str2double(temp{end}(end-(stubLngth-1+offset):end))+1);
        end
        temp{end}(end-(length(stub)-1):end) = stub;
        try
            StudyInstanceUID = strjoin(temp,'.');
        catch
            for l=1:numel(temp)
                if l==1
                    StudyInstanceUID = temp{l};
                else
                    StudyInstanceUID = [StudyInstanceUID '.' temp{l}];
                end
            end
        end
        meta.StudyInstanceUID = StudyInstanceUID;
        
        %transform slice pixel values to fit in UINT16
        intercept = 0;
        slope = 1;
        slice = vol(:,:,m,n);
        %negative values
        if min(slice(:))<0
            intercept = min(slice(:));
            slice = slice - intercept;
        end
        %floating point fractions
        temp = slice(:)*10;
        if any(mod(temp,10)~=0)
            slope = slope / 1e6;
            slice = slice * 1e6;
        end
        %large values
        if max(slice(:))>2^16-1
            slope = slope * max(slice(:)) / (2^16-1);
            slice = slice * (2^16-1) / max(slice(:));
        end
        slice = uint16(slice);
        
        meta.RescaleIntercept = intercept;
        meta.RescaleSlope = slope;
        
        %write with complete metadata, using 'Copy' argument
        dicomwrite(slice, [imageDir sliceName], meta, 'CreateMode', 'Copy')
        
        sliceCounter = sliceCounter + 1;
    end
end

if exist('h','var') && ishandle(h) && strcmp(get(h,'type'),'figure')
    close(h)
end

drawnow

end

function [S, imageFile] = parseInputs(varargin)

for n=1:nargin
    if ischar(varargin{n})
        imageFile = varargin{n};
    elseif isstruct(varargin{n})
        S = varargin{n};
    else
        error('Wrong input format.')
    end
end
if ~exist('imageFile','var')
    error('Must define image file name.')
end
if ~exist('S','var')
    error('Must input structure (see ''readImages'' function.')
end

end

function meta = getStaticDicomFields(meta, S)

%this function returns a structure of static dicom information

temp = S.studyDate;
if length(regexp(temp,':','split'))>1
    temp = regexp(temp,':','split');
    try
        temp = strjoin(temp,'');
    catch
        temp1 = temp;
        for n=1:numel(temp1)
            if n==1
                temp = temp1{n};
            else
                temp = [temp '' temp1{n}];
            end
        end
    end
end
meta.StudyDate = temp;
meta.SeriesDate = temp;
temp = S.studyTime;
if length(regexp(temp,':','split'))>1
    temp = regexp(temp,':','split');
    try
        temp = strjoin(temp,'');
    catch
        temp1 = temp;
        for n=1:numel(temp1)
            if n==1
                temp = temp1{n};
            else
                temp = [temp '' temp1{n}];
            end
        end
    end
end
meta.StudyTime = sprintf('%013.6f', str2double(temp));
meta.SeriesTime = sprintf('%013.6f', str2double(temp));
if isfield(S,'modality')
    meta.Modality = S.modality;
elseif isfield(S,'units')
    if strcmpi(S.units,'BQML')
        meta.Modality = 'PT';
    elseif strcmpi(S.units,'HU')
        meta.Modality = 'CT';
    end
end

if strcmpi(S.gantryModel,'1104')
    S.gantryModel = 'mCT';
elseif strcmpi(S.gantryModel,'2008')
    S.gantryModel = 'mMR';
end
meta.ManufacturerModelName = S.gantryModel;
if isfield(S,'patientWeight')
    meta.PatientWeight = S.patientWeight;
end
if isfield(S,'volume')
    meta.SliceThickness = abs(diff(S.bedRange3))/size(S.volume,3);
elseif isfield(S,'volumes')
    meta.SliceThickness = abs(diff(S.bedRange3))/size(S.volumes,3);
end
if isfield(S,'patientOrientation')
    meta.PatientPosition = S.patientOrientation;
end
meta.ImageOrientationPatient = S.IOP;
if isfield(S,'volume')
    meta.PixelSpacing = [abs(diff(S.bedRange2))/size(S.volume,2); ...
        abs(diff(S.bedRange1))/size(S.volume,1)];
elseif isfield(S,'volumes')
    meta.PixelSpacing = [abs(diff(S.bedRange2))/size(S.volumes,2); ...
        abs(diff(S.bedRange1))/size(S.volumes,1)];
end
if isfield(S,'tracer') %PET
    %tracer
    meta.RadiopharmaceuticalInformationSequence.Item_1.Radiopharmaceutical = S.tracer;
    meta.RadiopharmaceuticalInformationSequence.Item_1.RadionuclideTotalDose = S.injectedDose;
    temp = S.injectionDate;
    if length(regexp(temp,':','split'))>1
        temp = regexp(temp,':','split');
        try
            temp = strjoin(temp,'');
        catch
            temp1 = temp;
            for n=1:numel(temp1)
                if n==1
                    temp = temp1{n};
                else
                    temp = [temp '' temp1{n}];
                end
            end
        end
    end
    meta.RadiopharmaceuticalInformationSequence.Item_1.RadiopharmaceuticalStartDatetime = temp;
    temp = S.injectionTime;
    if length(regexp(temp,':','split'))>1
        temp = regexp(temp,':','split');
        try
            temp = strjoin(temp,'');
        catch
            temp1 = temp;
            for n=1:numel(temp1)
                if n==1
                    temp = temp1{n};
                else
                    temp = [temp '' temp1{n}];
                end
            end
        end
    end
    meta.RadiopharmaceuticalInformationSequence.Item_1.RadiopharmaceuticalStartTime = sprintf('%013.6f', str2double(temp));
    meta.StudyDescription = [meta.StudyDescription S.tracer];
end
if isfield(S,'volume')
    meta.NumberOfSlices = size(S.volume,3);
elseif isfield(S,'volumes')
    meta.NumberOfSlices = size(S.volumes,3);
end
if isfield(S,'units')
    meta.RescaleType = S.units;
    meta.Units = S.units;
end

meta.StudyDescription = [meta.StudyDescription '_' meta.Modality];
if isfield(S,'patientOrientation')
    meta.SeriesDescription = [meta.SeriesDescription '_' S.patientOrientation];
end
if isfield(S,'volume')
    meta.SeriesDescription = ...
        [meta.SeriesDescription '_' num2str(size(S.volume,1))];
elseif isfield(S,'volumes')
    meta.SeriesDescription = ...
        [meta.SeriesDescription '_' num2str(size(S.volumes,1))];
end

%generic
meta.SeriesNumber = 1;

end

function meta = getDynamicDicomFields(meta, S, n)

%this function modifies the dynamic dicom fields

if isfield(S,'endTimes')
    meta.ActualFrameDuration = (S.endTimes(n) - S.startTimes(n))*1e3;
    meta.FrameReferenceTime = mean([S.startTimes(n) S.endTimes(n)])*1e3; %maybe not accurate
end

%get acquistion time by adding frame start time to series time
srsTime = getSecondsFromTimeString(meta.SeriesTime);
if isfield(S,'startTimes')
    acqTime = srsTime + S.startTimes(n);
    if acqTime>=86400 %into the next day (never happens?)
        meta.AcquisitionDate = num2str(str2double(meta.SeriesDate) + 1);
        meta.AcquisitionTime = getTimeStringFromSeconds(acqTime-86400, 'dicom');
    else
        meta.AcquisitionDate = meta.SeriesDate;
        meta.AcquisitionTime = getTimeStringFromSeconds(acqTime, 'dicom');
    end
end

end

function out = initiateDicomStruct

out = struct;

out.Filename = '';
out.FileModDate = '';
out.FileSize = [];
out.Format = '';
out.FormatVersion = [];
out.Width = [];
out.Height = [];
out.BitDepth = [];
out.ColorType = '';
out.FileMetaInformationGroupLength = [];
out.FileMetaInformationVersion = [];
out.MediaStorageSOPClassUID = '';
out.MediaStorageSOPInstanceUID = '';
out.TransferSyntaxUID = '';
out.ImplementationClassUID = '';
out.ImplementationVersionName = '';
out.SOPClassUID = '';
out.SOPInstanceUID = '';
out.StudyDate = '';
out.SeriesDate = '';
out.AcquisitionDate = '';
out.ContentDate = '';
out.StudyTime = '';
out.SeriesTime = '';
out.AcquisitionTime = '';
out.ContentTime = '';
out.AccessionNumber = '';
out.Modality = '';
out.StudyDescription = '';
out.SeriesDescription = '';
out.ManufacturerModelName = '';
out.ConversionType = '';
out.ReferringPhysicianName = struct;
out.PatientName = struct;
out.PatientID = '';
out.PatientBirthDate = '';
out.PatientSex = '';
out.PatientWeight = [];
out.SliceThickness = [];
out.ActualFrameDuration = [];
out.PatientPosition = '';
out.SecondaryCaptureDeviceManufacturer = '';
out.SecondaryCaptureDeviceManufacturerModelName = '';
out.StudyInstanceUID = '';
out.SeriesInstanceUID = '';
out.StudyID = '';
out.SeriesNumber = [];
out.InstanceNumber = [];
out.ImagePositionPatient = [];
out.ImageOrientationPatient = [];
out.PatientOrientation = '';
out.SliceLocation = [];
out.SamplesPerPixel = [];
out.PhotometricInterpretation = '';
out.Rows = [];
out.Columns = [];
out.PixelSpacing = [];
out.BitsAllocated = [];
out.BitsStored = [];
out.HighBit = [];
out.PixelRepresentation = [];
out.SmallestImagePixelValue = [];
out.LargestImagePixelValue = [];
out.RescaleIntercept = [];
out.RescaleSlope = [];
out.RescaleType = '';
out.RadiopharmaceuticalInformationSequence = struct;
out.NumberOfSlices = [];
out.NumberOfTimeSlices = [];
out.Units = '';
out.FrameReferenceTime = [];

end

function out = getSOPClassUID(in)

if isfield(in,'modality')
    modality = in.modality;
else
    if isfield(in,'gantryModel')
        if strcmpi(in.gantryModel,'1104') || strcmpi(in.gantryModel(end-2:end),'mCT')
            modality = 'CT';
        elseif strcmpi(in.gantryModel,'2008') || strcmpi(in.gantryModel(end-2:end),'mMR')
            modality = 'MR';
        end
    end
    if isfield(in,'units') && strcmpi(in.units,'HU')
        modality = 'CT';
    end
    if isfield(in,'units') && strcmpi(in.units,'1/CM')
        modality = 'mu';
    end
    if isfield(in,'tracer') || ...
            (isfield(in,'units') && strcmpi(in.units,'BQML'))
        modality = 'PT';
    end
end

switch lower(modality)
    case 'pt'
        out = '1.2.840.10008.5.1.4.1.1.128';
    case 'mr'
        out = '1.2.840.10008.5.1.4.1.1.4';
    case 'ct'
        out = '1.2.840.10008.5.1.4.1.1.2';
    case 'mu'
        out = '1.2.840.10008.5.1.4.1.1.2';
    otherwise
        error(['SOP Class UID is not defined for ' modality '.'])
end

end

function UID = generateUID

%Matlab root to avoid conflicts with other vendors
ipt_root = '1.3.6.1.4.1.9590.100.1';

gUID = char(39);
for n = 1:39
    gUID(n) = num2str(randi(10)-1);
end

%remove leading zeros
while strcmp(gUID(1),'0')
    gUID(1) = [];
end

UID = [ipt_root '.2.' gUID];

end
