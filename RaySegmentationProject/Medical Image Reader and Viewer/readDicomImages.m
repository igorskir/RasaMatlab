function out = readDicomImages(in)

%this function accepts a single DICOM slice file
%and finds all slices in the directory which contribute to the same volume

%some dicom sub-rountines were adopted from
%  'DICOM to NIfTI converter, NIfTI tool and viewer' written by Xiangrui Li

%check core architecture for file path character
test = computer;
if strcmpi(test(1:2),'PC')
    pathChar = '\';
else
    pathChar = '/';
end

%vendor specific dictionary
dict = [];

if exist(in,'file')==7 %if passed a directory
    imageDir = in;
    if ~strcmpi(imageDir(end),pathChar), imageDir = [imageDir pathChar]; end
    %get all files in current directory
    test = dir([imageDir '*.*']);
    test([test.isdir]) = [];
    if isempty(test), error('No files in input directory'); end
elseif exist(in,'file')==2 %input is a specific file so extract path of containing folder
    temp = regexp(in,pathChar,'split');
    imageFile = temp{end};
    try
        imageDir = strjoin(temp(1:end-1),pathChar);
    catch
        for n=1:numel(temp)-1
            if n==1
                imageDir = temp{n};
            else
                imageDir = [imageDir pathChar temp{n}];
            end
        end
    end
    if isempty(imageDir), imageDir = cd; end
    if ~strcmpi(imageDir(end),pathChar), imageDir = [imageDir pathChar]; end
    if isempty(imageDir), imageDir = cd; end
    %get some volume attributes from slice
    [meta,~,dict] = dicm_hdr([imageDir pathChar imageFile]);
    if isempty(meta)
        disp('Input file is not recognized Dicom format')
        error('Input file is not recognized Dicom format')
    end
    if isfield(meta,'NumberOfSlices')
        nSlices = single(meta.NumberOfSlices);
    end
elseif exist(in,'file')==0
    error('Input must be file or directory')
end

%get all files in current directory
allFilesNOEXT = dir([imageDir '*.*']);
allFilesNOEXT([allFilesNOEXT.isdir]) = [];
allFiles = {allFilesNOEXT.name};

%wait bar display
startProg = 0;
if exist('nSlices','var')
    totalProg = nSlices;
else
    totalProg = length(allFiles);
end
h = waitbar(0,['Reading ' num2str(totalProg) ' files for Dicom data...']);

nSeries = 0;

oldNames = {};

if exist('imageFile','var') %single volume
    try
        if isfield(meta,'ImageType')
            imageType = meta.ImageType;
        end
        if isfield(meta,'ManufacturerModelName')
            gantryModel = meta.ManufacturerModelName;
        end
        if isfield(meta,'PatientPosition')
            patientOrientation = meta.PatientPosition;
        end
        if isfield(meta,'StudyDate')
            studyDate = meta.StudyDate;
        end
        if isfield(meta,'StudyTime')
            studyTime = meta.StudyTime;
        end
        if isfield(meta,'SeriesDate')
            seriesDate = meta.SeriesDate;
        end
        if isfield(meta,'SeriesTime')
            seriesTime = meta.SeriesTime;
        end
        if isfield(meta,'SeriesNumber')
            seriesNumber = single(meta.SeriesNumber);
        end
        if isfield(meta,'InstanceNumber')
            instanceNumber = single(meta.InstanceNumber);
        end
        if isfield(meta,'PatientWeight')
            patientWeight = single(meta.PatientWeight);
        end
        
        if isfield(meta,'NumberOfSlices')
            IPP = zeros(3,nSlices);
            volume = single(zeros(single(meta.Rows), single(meta.Columns), nSlices));
            instanceArr = zeros(nSlices,1)-1; %indexes the instance number
            sliceLoc = zeros(nSlices,1)-1e12; %tracks the slice location
            idx = zeros(nSlices,1); %indexes the file location
            %array of instance numbers to be included in the volume
            instances2include = floor((instanceNumber-1)/nSlices)*nSlices+(1:nSlices);
            %for CT, parameters are read for each slice, just initialize here
            if isfield(meta,'KVP')
                kVp = zeros(nSlices,1);
            end
            if isfield(meta,'XrayTubeCurrent')
                tubeCurrent = zeros(nSlices,1);
            end
            if isfield(meta,'ExposureTime')
                exposureTime = zeros(nSlices,1);
            end
            if isfield(meta,'SpiralPitchFactor')
                pitch = zeros(nSlices,1);
            end
            if isfield(meta,'Exposure')
                exposure = zeros(nSlices,1);
            end
        else
            IPP = zeros(3,1);
            volume = single(zeros(single(meta.Rows), single(meta.Columns)));
            instanceArr =  zeros(1,1)-1; %indexes the instance number
            sliceLoc = zeros(1,1)-1e12; %tracks the slice location
            idx = zeros(1,1); %indexes the file location
            %for CT, parameters are read for each slice, just initialize here
            if isfield(meta,'KVP')
                kVp = zeros(1,1);
            end
            if isfield(meta,'XrayTubeCurrent')
                tubeCurrent = zeros(1,1);
            end
            if isfield(meta,'ExposureTime')
                exposureTime = zeros(1,1);
            end
            if isfield(meta,'SpiralPitchFactor')
                pitch = zeros(1,1);
            end
            if isfield(meta,'Exposure')
                exposure = zeros(1,1);
            end
        end
        if isfield(meta,'Modality')
            modality = meta.Modality;
        end
        if isfield(meta,'RescaleType')
            units = meta.RescaleType;
        end
        if isfield(meta,'ImageOrientationPatient')
            IOP = meta.ImageOrientationPatient;
        end
        
        if isfield(meta,'RadiopharmaceuticalInformationSequence')
            if isfield(meta.RadiopharmaceuticalInformationSequence,'Item_1')
                if isfield(meta.RadiopharmaceuticalInformationSequence.Item_1,'Radiopharmaceutical')
                    tracer = meta.RadiopharmaceuticalInformationSequence.Item_1.Radiopharmaceutical;
                elseif isfield(meta.RadiopharmaceuticalInformationSequence.Item_1,'RadionuclideCodeSequence')
                    if isfield(meta.RadiopharmaceuticalInformationSequence.Item_1.RadionuclideCodeSequence,'Item_1')
                        if isfield(meta.RadiopharmaceuticalInformationSequence.Item_1.RadionuclideCodeSequence.Item_1,'CodeMeaning')
                            tracer = meta.RadiopharmaceuticalInformationSequence.Item_1.RadionuclideCodeSequence.Item_1.CodeMeaning;
                        end
                    end
                end
            end
            try
                injectionDate = meta.RadiopharmaceuticalInformationSequence.Item_1.RadiopharmaceuticalStartDatetime;
                injectionDate = injectionDate(1:8); %YYYYMMDD
            catch
                %                 warning('Injection date not found in PET DICOM data - using study date.')
                injectionDate = studyDate;
            end
            injectionTime = meta.RadiopharmaceuticalInformationSequence.Item_1.RadiopharmaceuticalStartTime;
            injectedDose = single(meta.RadiopharmaceuticalInformationSequence.Item_1.RadionuclideTotalDose);
            units = meta.Units;
            if isfield(meta, 'DoseCalibrationFactor')
                calibrationFctr = single(meta.DoseCalibrationFactor);
            end
            if isfield(meta.RadiopharmaceuticalInformationSequence.Item_1, 'RadionuclidePositronFraction')
                branchingFctr = single(meta.RadiopharmaceuticalInformationSequence.Item_1.RadionuclidePositronFraction);
            end
        end
        if isfield(meta,'NumberOfTimeSlices')
            nFrames = single(meta.NumberOfTimeSlices);
        end
        %this is not consistent for single multiple-bed volume (but it probably doesn't matter)
        if isfield(meta,'FrameReferenceTime') && isfield(meta,'ActualFrameDuration') %PET
            try
                startTime = getSecondsFromTimeString(meta.AcquisitionTime) - getSecondsFromTimeString(meta.SeriesTime);
            catch
                startTime = (single(meta.FrameReferenceTime)-single(meta.ActualFrameDuration)/2)/1000;
            end
            endTime = startTime + meta.ActualFrameDuration/1000;
            %                 %correct for delay between injection and acquisition times
            %                 if exist('injectionTime','var')
            %                     startTime = startTime + ...
            %                         (getSecondsFromTimeString(seriesTime) - getSecondsFromTimeString(injectionTime));
            %                     endTime = endTime + ...
            %                         (getSecondsFromTimeString(seriesTime) - getSecondsFromTimeString(injectionTime));
            %                 end
        end
        
        %this is to make the slice search efficient
        seed = find(strcmpi(imageFile,allFiles));
        %start with current file and move forward and backward from that point
        Z = seed; %current location
        ZR = seed; %this will track the reverse search progress
        ZF = seed; %this will track the forward search progress
        %get slice thickness
        if length(allFiles)>1 %it is best to use 'ImagePositionPatient' for slice spacing
            [temp1,~,dict] = dicm_hdr([imageDir allFiles{seed}],dict);
            if seed<length(allFiles)
                temp2 = dicm_hdr([imageDir allFiles{seed+1}],dict);
            else
                temp2 = [];
            end
            if isempty(temp2)
                temp2 = temp1;
                temp1 = dicm_hdr([imageDir allFiles{seed-1}],dict);
            end
            if ~isempty(temp1)
                %before calculating slice thickness, assure from same series
                check = strcmpi(temp1.ImageType,temp2.ImageType) && ...
                    strcmpi(temp1.SeriesInstanceUID,temp2.SeriesInstanceUID);
                if isfield(temp1,'StudyTime') && isfield(temp2,'StudyTime')
                    check = check && strcmpi(temp1.StudyTime,temp2.StudyTime);
                end
                if isfield(temp1,'SeriesNumber') && isfield(temp2,'SeriesNumber')
                    check = check && temp1.SeriesNumber==temp2.SeriesNumber;
                end
                if check
                    normVect = (temp2.ImagePositionPatient-temp1.ImagePositionPatient)/...
                        norm(temp2.ImagePositionPatient-temp1.ImagePositionPatient);
                end
            end
        end
        if ~exist('normVect','var')
            normVect = [0; 0; 1]; %arbitrary
        end
        %assure that adjacently parsed files increment only by slice thickness (for dynamic series)
        %             locR = normVect' * meta.ImagePositionPatient; %slice location boundary of reverse search
        %             locF = normVect' * meta.ImagePositionPatient; %slice location boundary of forward search
        currDir = 0; %this marks the current search direction (0-bkwrd, 1-frwrd)
        sliceCounter = 0;
        %set flag for long search warning
        warningFlag = 1;
        
        for n=1:length(allFiles)
            [temp,~,dict] = dicm_hdr([imageDir allFiles{Z}],dict);
            
            if ~isempty(temp) %dicom file
                sliceCond = true;
                
                %this is to select only one volume from dynamic PET series
                if exist('nSlices','var') && exist('nFrames','var')% && exist('ST','var')
                    sliceCond = ismember(temp.InstanceNumber,instances2include);
                    %                         %the following fails if incorrectly ordered slice files
                    %                         sliceCond = sliceCond & (round(abs(normVect'*temp.ImagePositionPatient-locR)*1e1)/1e1<=ST || ...
                    %                             round(abs(normVect'*temp.ImagePositionPatient-locF)*1e1)/1e1<=ST);
                end
                
                %condition for same volume and no repeat slices
                check = sliceCond && strcmpi(temp.ImageType,meta.ImageType) && ...
                    strcmpi(temp.SeriesInstanceUID,meta.SeriesInstanceUID);
                if check && isfield(temp,'StudyDate') && isfield(meta,'StudyDate')
                    check = strcmpi(temp.StudyDate,meta.StudyDate);
                end
                if check && isfield(temp,'StudyTime') && isfield(meta,'StudyTime')
                    check = strcmpi(temp.StudyTime,meta.StudyTime);
                end
                if check && isfield(temp,'SeriesNumber')
                    check = temp.SeriesNumber==seriesNumber;
                end
                if check && isfield(temp,'Rows')
                    check = size(volume,1)==single(temp.Rows);
                end
                if check && isfield(temp,'Columns')
                    check = size(volume,2)==single(temp.Columns);
                end
                if check && isfield(temp,'ImageOrientationPatient')
                    check = all(temp.ImageOrientationPatient==IOP);
                end
                if check && isfield(temp,'InstanceNumber')
                    check = ~ismember(temp.InstanceNumber,instanceArr);
                end
                if check && isfield(temp,'ImagePositionPatient')
                    check = ~ismember(normVect'*temp.ImagePositionPatient,sliceLoc);
                end
                
                if check
                    
                    sliceCounter = sliceCounter+1;
                    
                    idx(sliceCounter) = Z;
                    if isfield(temp,'InstanceNumber')
                        instanceArr(sliceCounter) = single(temp.InstanceNumber);
                    end
                    sliceLoc(sliceCounter) = normVect' * temp.ImagePositionPatient;
                    try
                        if isfield(temp,'RescaleSlope') && isfield(temp,'RescaleIntercept')
                            volume(:,:,sliceCounter) = single(dicm_img([imageDir allFiles{Z}]))*...
                                temp.RescaleSlope+temp.RescaleIntercept;
                        else
                            volume(:,:,sliceCounter) = single(dicm_img([imageDir allFiles{Z}]));
                        end
                    catch
                        error('Slices in same series may have inconsistent sizes')
                    end
                    %volume is ordered as (rows, columns, slices), i.e. X and Y are reversed
                    IPP(:,sliceCounter) = temp.ImagePositionPatient;
                    
                    %CT parameters may be different for each slice
                    if isfield(meta,'KVP')
                        kVp(sliceCounter) = single(meta.KVP);
                    end
                    if isfield(meta,'XrayTubeCurrent')
                        tubeCurrent(sliceCounter) = single(meta.XrayTubeCurrent);
                    end
                    if isfield(meta,'ExposureTime')
                        exposureTime(sliceCounter) = single(meta.ExposureTime);
                    end
                    if isfield(meta,'SpiralPitchFactor')
                        pitch(sliceCounter) = single(meta.SpiralPitchFactor);
                    end
                    if isfield(meta,'Exposure')
                        exposure(sliceCounter) = single(meta.Exposure);
                    end
                    
                    %smart search operations
                    if currDir
                        ZF = Z;
                        Z = ZF+1;
                        %                             locF = normVect' * temp.ImagePositionPatient;
                        %have you reached the last file?
                        if Z>length(allFiles)
                            currDir = 0;
                            Z = ZR-1;
                        end
                    else
                        ZR = Z;
                        Z = ZR-1;
                        %                             locR = normVect' * temp.ImagePositionPatient;
                        %have you reached the first file?
                        if Z<1
                            currDir = 1;
                            Z = ZF+1;
                        end
                    end
                else
                    %change search direction if file not found
                    if currDir
                        ZF = Z;
                        currDir = 0;
                        Z = ZR-1;
                        %have you reached the first file?
                        if Z<1
                            currDir = 1;
                            Z = ZF+1;
                        end
                    else
                        ZR = Z;
                        currDir = 1;
                        Z = ZF+1;
                        %have you reached the last file?
                        if Z>length(allFiles)
                            currDir = 0;
                            Z = ZR-1;
                        end
                    end
                end
                
                if exist('nSlices','var')
                    if exist('h','var') && ishandle(h) && strcmp(get(h,'type'),'figure')
                        waitbar(double((startProg+sliceCounter)/totalProg),h)
                    end
                    if sliceCounter==nSlices, break; end
                else
                    if exist('h','var') && ishandle(h) && strcmp(get(h,'type'),'figure')
                        waitbar(double((startProg+n)/totalProg),h)
                    end
                end
                %display text if taking longer than usual
                if exist('nSlices','var')
                    if n>4*nSlices && warningFlag
                        disp('Search is taking longer than usual, possibly missing files')
                        warningFlag = 0;
                    end
                end
            else  %non-dicom file, skip
                if currDir
                    ZF = Z;
                    Z = ZF+1;
                    %have you reached the last file?
                    if Z>length(allFiles)
                        currDir = 0;
                        Z = ZR-1;
                    end
                else
                    ZR = Z;
                    Z = ZR-1;
                    %have you reached the first file?
                    if Z<1
                        currDir = 1;
                        Z = ZF+1;
                    end
                end
            end
        end
        
        if exist('nSlices','var') && sliceCounter<nSlices
            if exist('h','var') && ishandle(h) && strcmp(get(h,'type'),'figure')
                close(h)
                drawnow
            end
            disp('Volume data are inconsistent or missing')
            error('Volume data are inconsistent or missing')
        end
        
        %use instance number instead to sort
        [test,idx] = sort(instanceArr);
        %             if ~all(diff(test)==1)
        %                 disp('Either duplicate data found or volume may be missing slices')
        %             end
        try
            volume(:,:,:) = volume(:,:,idx);
        catch
            disp('Slices in the current series may have different sizes')
            error('Slices in the current series may have different sizes')
        end
        IPP = IPP(:,idx);
        
        imSz1 = size(volume,1);
        imSz2 = size(volume,2);
        imSz3 = size(volume,3);
        
        %volume of X, Y, and Z coordinate values (patient space)
        %the following is equivalent to 'getVoxelLocs', but in DICOM space
        [temp2,temp1] = meshgrid(0:imSz2-1,0:imSz1-1);
        locX = zeros(imSz1,imSz2,imSz3);
        locY = zeros(imSz1,imSz2,imSz3);
        locZ = zeros(imSz1,imSz2,imSz3);
        temp = repmat(IPP(:,1),1,imSz1*imSz2)+...
            [IOP(1) IOP(2) IOP(3)]'*temp1(:)'*meta.PixelSpacing(2)+...
            [IOP(4) IOP(5) IOP(6)]'*temp2(:)'*meta.PixelSpacing(1);
        locX(:,:,1) = reshape(temp(1,:),imSz1,imSz2);
        locY(:,:,1) = reshape(temp(2,:),imSz1,imSz2);
        locZ(:,:,1) = reshape(temp(3,:),imSz1,imSz2);
        if imSz3>1
            temp = repmat(IPP(:,end),1,imSz1*imSz2)+...
                [IOP(1) IOP(2) IOP(3)]'*temp1(:)'*meta.PixelSpacing(2)+...
                [IOP(4) IOP(5) IOP(6)]'*temp2(:)'*meta.PixelSpacing(1);
            locX(:,:,end) = reshape(temp(1,:),imSz1,imSz2);
            locY(:,:,end) = reshape(temp(2,:),imSz1,imSz2);
            locZ(:,:,end) = reshape(temp(3,:),imSz1,imSz2);
            temp1 = [locX(1,1,1) locY(1,1,1) locZ(1,1,1)];
            temp2 = [locX(1,1,end) locY(1,1,end) locZ(1,1,end)];
            voxDim3 = round(norm(temp1-temp2)/(imSz3-1)*1e5)/1e5;
        else
            if isfield(meta,'SliceThickness')
                voxDim3 = meta.SliceThickness;
            else
                disp('Slice thickness cannot be determined')
            end
        end
        
        %get voxel dimensions in image (not necessarily scanner) space
        temp1 = [locX(1,1,1) locY(1,1,1) locZ(1,1,1)];
        temp2 = [locX(1,end,1) locY(1,end,1) locZ(1,end,1)];
        voxDim1 = round(norm(temp1-temp2)/(imSz2-1)*1e5)/1e5; %column spacing
        temp1 = [locX(1,1,1) locY(1,1,1) locZ(1,1,1)];
        temp2 = [locX(end,1,1) locY(end,1,1) locZ(end,1,1)];
        voxDim2 = round(norm(temp1-temp2)/(imSz1-1)*1e5)/1e5; %row spacing
        
        nSeries = nSeries+1;
        out.nSeries = nSeries;
        
        tempStruct = struct;
        
        if exist('seriesDate','var'), tempStruct.seriesDate = seriesDate; end
        if exist('seriesTime','var'), tempStruct.seriesTime = seriesTime; end
        if exist('modality','var'), tempStruct.modality = modality; end
        if exist('gantryModel','var'), tempStruct.gantryModel = gantryModel; end
        if exist('volume','var'), tempStruct.volumes = volume; end
        if exist('tracer','var'), tempStruct.tracer = tracer; end
        if exist('IOP','var'), tempStruct.IOP = IOP; end
        if exist('IPP','var'), tempStruct.IPP = IPP; end
        if exist('VPP','var'), tempStruct.VPP = VPP; end
        if exist('imSz1','var'), tempStruct.imSz1 = imSz1; end
        if exist('imSz2','var'), tempStruct.imSz2 = imSz2; end
        if exist('imSz3','var'), tempStruct.imSz3 = imSz3; end
        if exist('voxDim1','var'), tempStruct.voxDim1 = voxDim1; end
        if exist('voxDim2','var'), tempStruct.voxDim2 = voxDim2; end
        if exist('voxDim3','var'), tempStruct.voxDim3 = voxDim3; end
        if exist('injectedDose','var'), tempStruct.injectedDose = injectedDose; end
        if exist('patientOrientation','var'), tempStruct.patientOrientation = patientOrientation; end
        if exist('patientWeight','var'), tempStruct.patientWeight = patientWeight; end
        if exist('units','var'), tempStruct.units = units; end
        if exist('startTime','var'), tempStruct.startTimes = startTime; end
        if exist('endTime','var'), tempStruct.endTimes = endTime; end
        if exist('branchingFctr','var'), tempStruct.branchingFctr = branchingFctr; end
        if exist('calibrationFctr','var'), tempStruct.calibrationFctr = calibrationFctr; end
        %         if exist('studyDate','var'), tempStruct.studyDate = studyDate; end
        %         if exist('studyTime','var'), tempStruct.studyTime = studyTime; end
        if exist('injectionDate','var'), tempStruct.injectionDate = injectionDate; end
        if exist('injectionTime','var'), tempStruct.injectionTime = injectionTime; end
        
        if exist('kVp','var'), tempStruct.kVp = kVp(1); end
        if exist('tubeCurrent','var'), tempStruct.maxTubeCurrent = max(tubeCurrent); end
        if exist('tubeCurrent','var'), tempStruct.meanTubeCurrent = mean(tubeCurrent); end
        if exist('exposureTime','var'), tempStruct.exposureTime = exposureTime(1); end
        if exist('pitch','var'), tempStruct.pitch = pitch(1); end
        if exist('exposure','var'), tempStruct.maxExposure = max(exposure); end
        if exist('exposure','var'), tempStruct.meanExposure = mean(exposure); end
        
        % if exist('DICOMmeta','var'), tempStruct.DICOMmeta = DICOMmeta; end
        
        %save in under identifying series name
        if isfield(meta,'StudyDescription') && isfield(meta,'SeriesDescription')
            sName = [meta.StudyDescription '_' meta.SeriesDescription];
        else
            eval(['sName = [''Series_' num2str(nSeries) '''];'])
        end
        %remove illegal characters
        sName = rmIllChars(sName);
        out.(sName) = tempStruct;
    end
    
else %search all files and return every volume
    
    nFiles = length(allFiles);
    %initialize arrays to store identifying data
    imageTypeArr = cell(1,nFiles);
    studyDateArr = cell(1,nFiles);
    studyTimeArr = cell(1,nFiles);
    seriesDateArr = cell(1,nFiles);
    seriesTimeArr = cell(1,nFiles);
    seriesArr = zeros(1,nFiles)-1;
    seriesUIDArr = cell(1,nFiles);
    studyNameArr = cell(1,nFiles);
    seriesNameArr = cell(1,nFiles);
    instanceArr =  zeros(1,nFiles)-1; %indexes the instance number
    IPPArr = zeros(3,nFiles); %tracks the slice location
    sliceArr = cell(1,nFiles); %volume slice data
    rowsArr = zeros(1,nFiles);
    colsArr = zeros(1,nFiles);
    %optional values
    pixelSpacingArr = zeros(2,nFiles);
    sliceThicknessArr = zeros(1,nFiles);
    nSlicesArr = zeros(1,nFiles);
    modalityArr = cell(1,nFiles);
    gantryArr = cell(1,nFiles);
    unitsArr = cell(1,nFiles);
    patientOrientationArr = cell(1,nFiles);
    IOPArr = zeros(6,nFiles);
    patientWeightArr = zeros(1,nFiles);
    tracerArr = cell(1,nFiles);
    injDoseArr = zeros(1,nFiles);
    injDateArr = cell(1,nFiles);
    injTimeArr = cell(1,nFiles);
    calibrationFactorArr = zeros(1,nFiles);
    branchingFactorArr = zeros(1,nFiles);
    nFramesArr = zeros(1,nFiles);
    startTimesArr = zeros(1,nFiles)-1e12;
    endTimesArr = zeros(1,nFiles)-1e12;
    %CT
    kVpArr = zeros(1,nFiles);
    tubeCurrentArr = zeros(1,nFiles);
    exposureTimeArr = zeros(1,nFiles);
    pitchArr = zeros(1,nFiles);
    exposureArr = zeros(1,nFiles);
    %only used for spectroscopy
    VPPArr = zeros(3,nFiles);
    
    ignoreIdx = false(1,nFiles);
    
    for n=1:nFiles
        [meta,~,dict] = dicm_hdr([imageDir allFiles{n}],dict);
        
        if ~isempty(meta)
            if isfield(meta,'ImageType')
                imageTypeArr{n} = meta.ImageType;
            else
                imageTypeArr{n} = '';
            end
            if isfield(meta,'ManufacturerModelName')
                gantryArr{n} = meta.ManufacturerModelName;
            end
            if isfield(meta,'PatientPosition')
                patientOrientationArr{n} = meta.PatientPosition;
            end
            if isfield(meta,'SeriesDate')
                seriesDateArr{n} = meta.SeriesDate;
            end
            if isfield(meta,'SeriesTime')
                seriesTimeArr{n} = meta.SeriesTime;
            end
            if isfield(meta,'SeriesNumber')
                seriesArr(n) = meta.SeriesNumber;
            end
            if isfield(meta,'StudyDate')
                studyDateArr{n} = meta.StudyDate;
            end
            if isfield(meta,'StudyTime')
                studyTimeArr{n} = meta.StudyTime;
            end
            if isfield(meta,'SeriesInstanceUID')
                seriesUIDArr{n} = meta.SeriesInstanceUID;
            end
            if isfield(meta,'InstanceNumber')
                instanceArr(n) = meta.InstanceNumber;
            end
            if isfield(meta,'StudyDescription')
                studyNameArr{n} = meta.StudyDescription;
            end
            if isfield(meta,'SeriesDescription')
                seriesNameArr{n} = meta.SeriesDescription;
            end
            if isfield(meta,'PatientWeight')
                patientWeightArr(n) = single(meta.PatientWeight);
            end
            
            %process some image and spectroscopy data differently
            if exist('spectro','var')
                rowsArr(n) = spectro.Rows;
                colsArr(n) = spectro.Columns;
                nSlicesArr(n) = spectro.NumberOfFrames;
                if colsArr(n)==1 && rowsArr(n)==1 && nSlicesArr(n)==1
                    modalityArr{n} = 'SVS';
                else
                    modalityArr{n} = 'CSI';
                end
                pixelSpacingArr(:,n) = single(spectro.PixelSpacing);
                sliceThicknessArr(n) = single(spectro.SliceThickness);
                IOPArr(:,n) =  spectro.ImageOrientationPatient;
                IPPArr(:,n) = spectro.ImagePositionPatient;
                VPPArr(:,n) = spectro.VoiPosition;
            else
                if isfield(meta,'RescaleSlope') && isfield(meta,'RescaleIntercept')
                    sliceArr{n} = single(dicm_img([imageDir allFiles{n}]))*...
                        meta.RescaleSlope+meta.RescaleIntercept;
                else
                    sliceArr{n} = dicm_img([imageDir allFiles{n}]);
                end
                if isfield(meta,'ImagePositionPatient')
                    IPPArr(:,n) = meta.ImagePositionPatient;
                end
                if isfield(meta,'ImageOrientationPatient')
                    IOPArr(:,n) = meta.ImageOrientationPatient;
                end
                %volume is ordered as (rows, columns, slices), i.e. X and Y are reversed
                rowsArr(n) = meta.Rows;
                colsArr(n) = meta.Columns;
                
                if isfield(meta,'PixelSpacing')
                    pixelSpacingArr(:,n) = single(meta.PixelSpacing);
                end
                if isfield(meta,'SliceThickness')
                    sliceThicknessArr(n) = single(meta.SliceThickness);
                end
                if isfield(meta,'NumberOfSlices')
                    nSlicesArr(n) = meta.NumberOfSlices;
                end
                if isfield(meta,'Modality')
                    modalityArr{n} = meta.Modality;
                end
                if isfield(meta,'RescaleType')
                    unitsArr{n} = meta.RescaleType;
                end
                if isfield(meta,'RadiopharmaceuticalInformationSequence')
                    if isfield(meta.RadiopharmaceuticalInformationSequence,'Item_1')
                        if isfield(meta.RadiopharmaceuticalInformationSequence.Item_1,'Radiopharmaceutical')
                            tracerArr{n} = meta.RadiopharmaceuticalInformationSequence.Item_1.Radiopharmaceutical;
                        elseif isfield(meta.RadiopharmaceuticalInformationSequence.Item_1,'RadionuclideCodeSequence')
                            if isfield(meta.RadiopharmaceuticalInformationSequence.Item_1.RadionuclideCodeSequence,'Item_1')
                                if isfield(meta.RadiopharmaceuticalInformationSequence.Item_1.RadionuclideCodeSequence.Item_1,'CodeMeaning')
                                    tracerArr{n} = meta.RadiopharmaceuticalInformationSequence.Item_1.RadionuclideCodeSequence.Item_1.CodeMeaning;
                                end
                            end
                        end
                    end
                    try
                        injDateArr{n} = meta.RadiopharmaceuticalInformationSequence.Item_1.RadiopharmaceuticalStartDatetime;
                        injDateArr{n} = injDateArr{n}(1:8); %YYYYMMDD
                    catch
                        %                     warning('Injection date not found in DICOM data. Using study date.')
                        injDateArr{n} = studyDateArr{n};
                    end
                    injTimeArr{n} = meta.RadiopharmaceuticalInformationSequence.Item_1.RadiopharmaceuticalStartTime;
                    injDoseArr(n) = single(meta.RadiopharmaceuticalInformationSequence.Item_1.RadionuclideTotalDose);
                    unitsArr{n} = meta.Units;
                    if isfield(meta, 'DoseCalibrationFactor')
                        calibrationFactorArr(n) = single(meta.DoseCalibrationFactor);
                    end
                    if isfield(meta.RadiopharmaceuticalInformationSequence.Item_1, 'RadionuclidePositronFraction')
                        branchingFactorArr(n) = single(meta.RadiopharmaceuticalInformationSequence.Item_1.RadionuclidePositronFraction);
                    end
                end
                if isfield(meta,'NumberOfTimeSlices')
                    nFramesArr(n) = meta.NumberOfTimeSlices;
                end
                %this is not consistent for single multiple-bed volume (but it probably doesn't matter)
                if isfield(meta,'FrameReferenceTime') && isfield(meta,'ActualFrameDuration') %PET
                    try
                        startTimesArr(n) = getSecondsFromTimeString(meta.AcquisitionTime) - getSecondsFromTimeString(meta.SeriesTime);
                    catch
                        startTimesArr(n) = (meta.FrameReferenceTime - meta.ActualFrameDuration/2)/1000;
                    end
                    endTimesArr(n) = startTimesArr(n) + meta.ActualFrameDuration/1000;
                    %                     %correct for delay between injection and acquisition times
                    %                     if ~isempty(injTimeArr{n})
                    %                         startTimesArr(n) = startTimesArr(n) + ...
                    %                             (getSecondsFromTimeString(seriesTimeArr{n}) - getSecondsFromTimeString(injTimeArr{n}));
                    %                         endTimesArr(n) = endTimesArr(n) + ...
                    %                             (getSecondsFromTimeString(seriesTimeArr{n}) - getSecondsFromTimeString(injTimeArr{n}));
                    %                     end
                end
                %for CT
                if isfield(meta,'KVP')
                    kVpArr(n) = single(meta.KVP);
                end
                if isfield(meta,'XrayTubeCurrent')
                    tubeCurrentArr(n) = single(meta.XrayTubeCurrent);
                end
                if isfield(meta,'ExposureTime')
                    exposureTimeArr(n) = single(meta.ExposureTime);
                end
                if isfield(meta,'SpiralPitchFactor')
                    pitchArr(n) = single(meta.SpiralPitchFactor);
                end
                if isfield(meta,'Exposure')
                    exposureArr(n) = single(meta.Exposure);
                end
            end
        else
            %             disp(['''' [imageDir allFiles{n}] ''' is not Dicom format and will be skipped'])
            ignoreIdx(n) = true;
        end
        if exist('h','var') && ishandle(h) && strcmp(get(h,'type'),'figure')
            waitbar(double((startProg+n)/totalProg),h)
        end
    end
    
    %remove ignored indeces
    imageTypeArr(ignoreIdx) = [];
    studyDateArr(ignoreIdx) = [];
    studyTimeArr(ignoreIdx) = [];
    seriesDateArr(ignoreIdx) = [];
    seriesTimeArr(ignoreIdx) = [];
    seriesArr(ignoreIdx) = [];
    studyNameArr(ignoreIdx) = [];
    seriesNameArr(ignoreIdx) = [];
    seriesUIDArr(ignoreIdx) = [];
    instanceArr(ignoreIdx) = [];
    IPPArr(:,ignoreIdx) = [];
    IOPArr(:,ignoreIdx) = [];
    sliceArr(ignoreIdx) = [];
    rowsArr(ignoreIdx) = [];
    colsArr(ignoreIdx) = [];
    pixelSpacingArr(:,ignoreIdx) = [];
    sliceThicknessArr(ignoreIdx) = [];
    nSlicesArr(ignoreIdx) = [];
    modalityArr(ignoreIdx) = [];
    gantryArr(ignoreIdx) = [];
    unitsArr(ignoreIdx) = [];
    patientOrientationArr(ignoreIdx) = [];
    patientWeightArr(ignoreIdx) = [];
    tracerArr(ignoreIdx) = [];
    injDoseArr(ignoreIdx) = [];
    injDateArr(ignoreIdx) = [];
    injTimeArr(ignoreIdx) = [];
    calibrationFactorArr(ignoreIdx) = [];
    branchingFactorArr(ignoreIdx) = [];
    nFramesArr(ignoreIdx) = [];
    startTimesArr(ignoreIdx) = [];
    endTimesArr(ignoreIdx) = [];
    VPPArr(:,ignoreIdx) = [];
    kVpArr(ignoreIdx) = [];
    tubeCurrentArr(ignoreIdx) = [];
    exposureTimeArr(ignoreIdx) = [];
    pitchArr(ignoreIdx) = [];
    exposureArr(ignoreIdx) = [];
    
    %identify all series, leave the last 2 columns of 'seriesData' if adding more identifying fields
    seriesData = cell(1,2); %series name and file indeces
    seriesUIDs = unique(seriesUIDArr);
    for n=1:numel(seriesUIDs)
        idx = strcmpi(seriesUIDArr,seriesUIDs(n));
        imageTypes = unique(imageTypeArr(idx));
        for m=1:numel(imageTypes)
            idx = strcmpi(seriesUIDArr,seriesUIDs(n)) & ...
                strcmpi(imageTypeArr,imageTypes(m));
            unique(pixelSpacingArr(:,idx)','rows');
            rowscols = unique([rowsArr(idx); colsArr(idx)]','rows');
            for l=1:numel(rowscols)/2
                idx = strcmpi(seriesUIDArr,seriesUIDs(n)) & ...
                    strcmpi(imageTypeArr,imageTypes(m)) & ...
                    (rowsArr==rowscols(l,1) & colsArr==rowscols(l,2));
                pixelSpacings = unique(pixelSpacingArr(:,idx)','rows');
                for k=1:numel(pixelSpacings)/2
                    idx = strcmpi(seriesUIDArr,seriesUIDs(n)) & ...
                        strcmpi(imageTypeArr,imageTypes(m)) & ...
                        (rowsArr==rowscols(l,1) & colsArr==rowscols(l,2)) & ...
                        (pixelSpacingArr(1,:)==pixelSpacings(k,1) & ...
                        pixelSpacingArr(2,:)==pixelSpacings(k,2));
                    IOPs = unique(IOPArr(:,idx)','rows');
                    for j=1:numel(IOPs)/6
                        idx = strcmpi(seriesUIDArr,seriesUIDs(n)) & ...
                            strcmpi(imageTypeArr,imageTypes(m)) & ...
                            (rowsArr==rowscols(l,1) & colsArr==rowscols(l,2)) & ...
                            (pixelSpacingArr(1,:)==pixelSpacings(k,1) & ...
                            pixelSpacingArr(2,:)==pixelSpacings(k,2)) & ...
                            (IOPArr(1,:)==IOPs(j,1) & ...
                            IOPArr(2,:)==IOPs(j,2) & ...
                            IOPArr(3,:)==IOPs(j,3) & ...
                            IOPArr(4,:)==IOPs(j,4) & ...
                            IOPArr(5,:)==IOPs(j,5) & ...
                            IOPArr(6,:)==IOPs(j,6));
                        
                        idx = find(idx);
                        nSeries = nSeries+1;
                        
                        if ~isempty(studyNameArr(idx(1))) || ~isempty(seriesNameArr(idx(1)))
                            if ~isempty(studyNameArr(idx(1))) && ~isempty(seriesNameArr(idx(1)))
                                seriesName = {[studyNameArr{idx(1)} '_' seriesNameArr{idx(1)}]};
                            else
                                if ~isempty(studyNameArr(idx(1)))
                                    seriesName = {studyNameArr{idx(1)}};
                                else
                                    seriesName = {seriesNameArr{idx(1)}};
                                end
                            end
                        else
                            eval(['seriesName = {''Series_' num2str(nSeries) '''};'])
                        end
                        %sort according to instance number (IMPORTANT)
                        [~,srt] = sort(instanceArr(idx));
                        idx = idx(srt);
                        %save data in cell matrix
                        seriesData(nSeries,:) = [seriesName {idx}];
                    end
                end
            end
        end
    end
    
    out.nSeries = nSeries;
    
    %now save all volumes in output
    for n=1:nSeries
        try %assure valid and complete series
            %identify dublicate series if any, index data are sequentially sorted by instance number
            idx = seriesData{n,end};
            testIdx = [seriesData{n,end}(diff(instanceArr(seriesData{n,end}))>0) seriesData{n,end}(end)];
            if length(idx)>1
                if length(testIdx)<length(idx)
                    disp(['Duplicate data found for series: ' seriesData{n,end-1}])
                    idx = testIdx;
                    seriesData{n,end} = idx;
                end
            end
            %get series variables
            if ~isempty(modalityArr{idx(1)})
                modality = modalityArr{idx(1)};
            else
                clear modality
            end
            if ~isempty(gantryArr{idx(1)})
                gantryModel = gantryArr{idx(1)};
            else
                clear gantryModel
            end
            if ~isempty(tracerArr{idx(1)})
                tracer = tracerArr{idx(1)};
            else
                clear tracer
            end
            if injDoseArr(idx(1))>0
                injectedDose = injDoseArr(idx(1));
            else
                clear injectedDose
            end
            if ~isempty(patientOrientationArr{idx(1)})
                patientOrientation = patientOrientationArr{idx(1)};
            else
                clear patientOrientation
            end
            if patientWeightArr(idx(1))>0
                patientWeight = patientWeightArr(idx(1));
            else
                clear patientWeight
            end
            if ~isempty(unitsArr{idx(1)})
                units = unitsArr{idx(1)};
            else
                clear units
            end
            if branchingFactorArr(idx(1))>0
                branchingFctr = branchingFactorArr(idx(1));
            else
                clear branchingFctr
            end
            if calibrationFactorArr(idx(1))>0
                calibrationFctr = calibrationFactorArr(idx(1));
            else
                clear calibrationFctr
            end
            if ~isempty(injDateArr{idx(1)})
                injectionDate = injDateArr{idx(1)};
            else
                clear injectionDate
            end
            if ~isempty(injTimeArr{idx(1)})
                injectionTime = injTimeArr{idx(1)};
            else
                clear injectionTime
            end
            %CT parameters may change for each slice
            if kVpArr(idx(1))>0
                kVp = mean(kVpArr(idx(:))); %should be constant
            else
                clear kVp
            end
            if tubeCurrentArr(idx(1))>0
                maxTubeCurrent = max(tubeCurrentArr(idx(:)));
                meanTubeCurrent = mean(tubeCurrentArr(idx(:)));
            else
                clear maxTubeCurrent meanTubeCurrent
            end
            if exposureTimeArr(idx(1))>0
                exposureTime = mean(exposureTimeArr(idx(:))); %should be constant
            else
                clear exposureTime
            end
            if pitchArr(idx(1))>0
                pitch = mean(pitchArr(idx(:))); %should be constant
            else
                clear pitch
            end
            if exposureArr(idx(1))>0
                maxExposure = max(exposureArr(idx(:)));
                meanExposure = mean(exposureArr(idx(:)));
            else
                clear maxExposure meanExposure
            end
            
            %build IPP array
            IPP = IPPArr(:,idx);
            
            if ~strcmpi(modality,'SVS') && ~strcmpi(modality,'CSI')
                %check for duplicate slice locations in same series, assume time frames
                unq1 = unique(IPP(1,:));
                unq2 = unique(IPP(2,:));
                unq3 = unique(IPP(3,:));
                nUnq = max([length(unq1) length(unq2) length(unq3)]);
                if nUnq<length(idx)
                    %multiple time frames found
                    if mod(length(idx),nUnq)~=0
                        %incomplete series?
                        if exist('h','var') && ishandle(h) && strcmp(get(h,'type'),'figure')
                            close(h)
                            drawnow
                        end
                        disp(['Incomplete dynamic series: ' seriesData{n,end-1}])
                        error(['Incomplete dynamic series: ' seriesData{n,end-1}])
                    end
                    nFrames = length(idx)/nUnq;
                else
                    if nSlicesArr(idx(1))>0 && length(idx)<nSlicesArr(idx(1))
                        if exist('h','var') && ishandle(h) && strcmp(get(h,'type'),'figure')
                            close(h)
                            drawnow
                        end
                        disp(['Volume is missing slices: ' seriesData{n,end-1}])
                        error(['Volume is missing slices: ' seriesData{n,end-1}])
                    end
                    nFrames = 1;
                end
            end
            
            %volume sizes
            imSz1 = rowsArr(idx(1));
            imSz2 = colsArr(idx(1));
            if nSlicesArr(idx(1))>0
                imSz3 = nSlicesArr(idx(1));
            else
                imSz3 = length(idx)/nFrames;
            end
            
            %frame times
            if any(startTimesArr(idx)~=-1e12)
                startTimes = zeros(1,nFrames);
                endTimes = zeros(1,nFrames);
                for m=1:nFrames
                    startTimes(m) = startTimesArr((m-1)*imSz3+1);
                    endTimes(m) = endTimesArr((m-1)*imSz3+1);
                end
            else
                clear startTimes endTimes
            end
            
            %adjust IPP
            if ~strcmpi(modality,'SVS') && ~strcmpi(modality,'CSI')
                IPP = IPP(:,1:imSz3);
            end
            %IOP
            IOP = IOPArr(:,idx(1));
            
            %build volume (not for spectroscopy data)
            if ~isempty(sliceArr{idx(1)})
                try
                    volumes = single(reshape([sliceArr{idx}],imSz1,imSz2,imSz3,nFrames));
                catch
                    disp('Slices belonging to one series in the current directory may have different sizes')
                    error('Slices belonging to one series in the current directory are inconsistent')
                end
                
                %volume of X, Y, and Z coordinate values (patient space)
                %the following is equivalent to 'getVoxelLocs', but in DICOM space
                [temp2,temp1] = meshgrid(0:imSz2-1,0:imSz1-1);
                locX = zeros(imSz1,imSz2,imSz3);
                locY = zeros(imSz1,imSz2,imSz3);
                locZ = zeros(imSz1,imSz2,imSz3);
                ps1 = pixelSpacingArr(1,idx(1));
                ps2 = pixelSpacingArr(2,idx(1));
                st = sliceThicknessArr(idx(1));
                temp = repmat(IPP(:,1),1,imSz1*imSz2)+...
                    [IOP(1) IOP(2) IOP(3)]'*temp1(:)'*ps1+...
                    [IOP(4) IOP(5) IOP(6)]'*temp2(:)'*ps2;
                locX(:,:,1) = reshape(temp(1,:),imSz1,imSz2);
                locY(:,:,1) = reshape(temp(2,:),imSz1,imSz2);
                locZ(:,:,1) = reshape(temp(3,:),imSz1,imSz2);
                voxDim3 = st;
                if imSz3>1
                    temp = repmat(IPP(:,end),1,imSz1*imSz2)+...
                        [IOP(1) IOP(2) IOP(3)]'*temp1(:)'*ps1+...
                        [IOP(4) IOP(5) IOP(6)]'*temp2(:)'*ps2;
                    locX(:,:,end) = reshape(temp(1,:),imSz1,imSz2);
                    locY(:,:,end) = reshape(temp(2,:),imSz1,imSz2);
                    locZ(:,:,end) = reshape(temp(3,:),imSz1,imSz2);
                    temp1 = [locX(1,1,1) locY(1,1,1) locZ(1,1,1)];
                    temp2 = [locX(1,1,end) locY(1,1,end) locZ(1,1,end)];
                    voxDim3 = round(norm(temp1-temp2)/(imSz3-1)*1e5)/1e5;
                end
                
                %get voxel dimensions in image (not necessarily scanner) space
                temp1 = [locX(1,1,1) locY(1,1,1) locZ(1,1,1)];
                temp2 = [locX(1,end,1) locY(1,end,1) locZ(1,end,1)];
                voxDim1 = round(norm(temp1-temp2)/(imSz2-1)*1e5)/1e5; %column spacing
                temp1 = [locX(1,1,1) locY(1,1,1) locZ(1,1,1)];
                temp2 = [locX(end,1,1) locY(end,1,1) locZ(end,1,1)];
                voxDim2 = round(norm(temp1-temp2)/(imSz1-1)*1e5)/1e5; %row spacing
            else
                %spectroscopy data (always one 'slice')
                voxDim1 =  pixelSpacingArr(2,idx(1));
                voxDim2 =  pixelSpacingArr(1,idx(1));
                voxDim3 =  sliceThicknessArr(idx(1));
                VPP = VPPArr(:,idx);
            end
            
            %initialize structure
            S = struct;
            
            if isfield(meta,'SeriesDate')
                S.seriesDate = meta.SeriesDate;
            elseif isfield(meta,'StudyDate')
                S.seriesDate = meta.StudyDate;
            end
            if isfield(meta,'SeriesTime')
                S.seriesTime = meta.SeriesTime;
            elseif isfield(meta,'StudyTime')
                S.seriesTime = meta.StudyTime;
            end
            
            if exist('modality','var'), S.modality = modality; end
            if exist('gantryModel','var'), S.gantryModel = gantryModel; end
            if exist('volumes','var'), S.volumes = volumes; end
            if exist('tracer','var'), S.tracer = tracer; end
            if exist('IOP','var'), S.IOP = IOP; end
            if exist('IPP','var'), S.IPP = IPP; end
            if exist('VPP','var'), S.VPP = VPP; end
            if exist('imSz1','var'), S.imSz1 = imSz1; end
            if exist('imSz2','var'), S.imSz2 = imSz2; end
            if exist('imSz3','var'), S.imSz3 = imSz3; end
            if exist('voxDim1','var'), S.voxDim1 = voxDim1; end
            if exist('voxDim2','var'), S.voxDim2 = voxDim2; end
            if exist('voxDim3','var'), S.voxDim3 = voxDim3; end
            if exist('injectedDose','var'), S.injectedDose = injectedDose; end
            if exist('patientOrientation','var'), S.patientOrientation = patientOrientation; end
            if exist('patientWeight','var'), S.patientWeight = patientWeight; end
            if exist('units','var'), S.units = units; end
            if exist('startTimes','var'), S.startTimes = startTimes; end
            if exist('endTimes','var'), S.endTimes = endTimes; end
            if exist('branchingFctr','var'), S.branchingFctr = branchingFctr; end
            if exist('calibrationFctr','var'), S.calibrationFctr = calibrationFctr; end
            %PET
            if exist('injectionDate','var'), S.injectionDate = injectionDate; end
            if exist('injectionTime','var'), S.injectionTime = injectionTime; end
            %CT
            if exist('kVp','var'), S.kVp = kVp; end
            if exist('maxTubeCurrent','var'), S.maxTubeCurrent = maxTubeCurrent; end
            if exist('meanTubeCurrent','var'), S.meanTubeCurrent = meanTubeCurrent; end
            if exist('exposureTime','var'), S.exposureTime = exposureTime; end
            if exist('pitch','var'), S.pitch = pitch; end
            if exist('maxExposure','var'), S.maxExposure = maxExposure; end
            if exist('meanExposure','var'), S.meanExposure = meanExposure; end
            
            sName = seriesData{n,end-1};
            
            %remove illegal characters
            sName = rmIllChars(sName);
            
            if isfield(S,'volumes') && size(S.volumes,4)>1
                sName = [sName '_' num2str(size(S.volumes,4)) 'frames'];
            end
            %correct for unnamed series
            if strcmpi(sName,'_')
                eval(['sName = [''Series_' num2str(n) '''];'])
            end
            %character limit
            if numel(sName)>60, sName = sName(1:60); end
            %check if series name already exists in output structure
            temp = sName;
            fn = fieldnames(out);
            cnt = 2; %it's 2 because first will be numbered later
            while ismember(temp,fn)
                %set number to old series
                oldNames = [oldNames {sName}];
                temp = [sName '_' num2str(cnt)];
                cnt = cnt+1;
            end
            sName = temp;
            if strcmpi(sName(1),'_'), sName = sName(2:end); end
            
            out.(sName) = S;
            
        catch
            out.nSeries = out.nSeries-1;
        end
    end
end

%if repeated series found, set number 1 to name of first
oldNames = unique(oldNames);
fn = fieldnames(out);
v = struct2cell(out);
for n=1:length(oldNames)
    eval(['fn{strcmpi(''' oldNames{n} ''',fn)} = ''' oldNames{n} '_1'';'])
    out = cell2struct(v,fn);
end

if exist('h','var') && ishandle(h) && strcmp(get(h,'type'),'figure')
    close(h)
    drawnow
end

end

function sName = rmIllChars(sName)

try
    sName = strjoin(regexp(sName,'[!@#$%^&*+./''";:<>=|]','split'),'_');
    sName = strjoin(regexp(sName,'\','split'),'_');
    sName = strjoin(regexp(sName,'/','split'),'_');
    sName = strjoin(regexp(sName,'-','split'),'_');
    sName = strjoin(regexp(sName,'(','split'),'_');
    sName = strjoin(regexp(sName,')','split'),'_');
    sName = strjoin(regexp(sName,'[','split'),'_');
    sName = strjoin(regexp(sName,']','split'),'_');
    sName = strjoin(regexp(sName,',','split'),'_');
    sName = strjoin(regexp(sName,'`','split'),'_');
    sName = strjoin(regexp(sName,' ','split'),'_');
catch
    temp = regexp(sName,'[!@#$%^&*+./''";:<>=|]','split');
    for n=1:numel(temp)
        if n==1
            sName = temp{n};
        else
            sName = [sName '_' temp{n}];
        end
    end
    temp = regexp(sName,'\','split');
    for n=1:numel(temp)
        if n==1
            sName = temp{n};
        else
            sName = [sName '_' temp{n}];
        end
    end
    temp = regexp(sName,'/','split');
    for n=1:numel(temp)
        if n==1
            sName = temp{n};
        else
            sName = [sName '_' temp{n}];
        end
    end
    temp = regexp(sName,'-','split');
    for n=1:numel(temp)
        if n==1
            sName = temp{n};
        else
            sName = [sName '_' temp{n}];
        end
    end
    temp = regexp(sName,'(','split');
    for n=1:numel(temp)
        if n==1
            sName = temp{n};
        else
            sName = [sName '_' temp{n}];
        end
    end
    temp = regexp(sName,')','split');
    for n=1:numel(temp)
        if n==1
            sName = temp{n};
        else
            sName = [sName '_' temp{n}];
        end
    end
    temp = regexp(sName,'[','split');
    for n=1:numel(temp)
        if n==1
            sName = temp{n};
        else
            sName = [sName '_' temp{n}];
        end
    end
    temp = regexp(sName,']','split');
    for n=1:numel(temp)
        if n==1
            sName = temp{n};
        else
            sName = [sName '_' temp{n}];
        end
    end
    temp = regexp(sName,',','split');
    for n=1:numel(temp)
        if n==1
            sName = temp{n};
        else
            sName = [sName '_' temp{n}];
        end
    end
    temp = regexp(sName,'`','split');
    for n=1:numel(temp)
        if n==1
            sName = temp{n};
        else
            sName = [sName '_' temp{n}];
        end
    end
    temp = regexp(sName,' ','split');
    for n=1:numel(temp)
        if n==1
            sName = temp{n};
        else
            sName = [sName '_' temp{n}];
        end
    end
end

%remove redundant underscores
test = regexp(sName,'_','split');
sName = strjoin(test(~cellfun(@isempty, test)),'_');

end