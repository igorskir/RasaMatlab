
%% Import dicom images
% The function 'readImages' reads dicom image data from an image file or
% folder. Important attributes are stored in a convenient structure, which
% is used as the input for many other MATLAB functions. It can be called
% with no arguments, which will prompt for file selection, or a specific
% image file or folder. Dynamic series are loaded by passing only the
% series directory folder, with no file name.
%
% Output volume X and Y axes count columns and rows, respectively. This is
% intuitive, but is opposite dicom convention
%
% NOTE: Output bed ranges are in image coordinates which are not scanner
% coordinates (unless IOP==[1; 0; 0; 0; 1; 0])
%
%% Viewing single images
%
% Read PET data from dicom directory
PET = readImages('C:\PET\');
% PET volume is already aligned to scanner coordinate axes
tic; VolumeViewer3D(PET); toc
%% Viewing fused images
%
% Read CT data from dicom directory
CT = readImages('C:\CT\');
% Both volumes are interpolated to be displayed on a common voxel grid
% The first input image is the target orientation, so the following call
% rotates PET into CT orientation, which is anti-aligned along the Z axis
tic; VolumeViewer3D(CT,PET); toc
% Pass flag to align both volumes to scanner coordinate axes
% In this case, this is equivalent to calling VolumeViewer3D(PET,CT)
tic; VolumeViewer3D(CT,PET,'align'); toc
% Can pass flag to show only shared space between both volumes
% Smaller matrices result in faster interpolation
tic; VolumeViewer3D(PET,CT,'trim'); toc
% Similar for PET/MR, where the relative orientations may be very different
MR = readImages('C:\MR\');
PET_MR = readImages('C:\PET(MR)\');
VolumeViewer3D(MR,PET_MR)
%% Drawing ROI masks
%
% Region-of-interest masks can be generated with the argument 'drawROI'
% All other arguments are ignored, and if a second volume is passed, it
% will be defined in the space of the first, as will the output mask
mask = VolumeViewer3D(PET,CT,'drawROI');
%% Output dicom images
%
% Write dicom PET slices to new directory
writeDicomImage(PET,'C:\newPET\')