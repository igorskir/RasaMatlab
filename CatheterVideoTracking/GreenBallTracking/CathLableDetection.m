clc
close all

%% Open video file
[filename, pathname] = ...
     uigetfile({'*.mp4';'*.avi';'*.mpeg'},'File selector');

%% Reader and player configuration
videoReader = VideoReader(fullfile(pathname, filename));
videoPlayer = vision.VideoPlayer('Position', [20, 400, 700, 400]);

%% Get number of frames from video
nFrames = round(videoReader.Duration * videoReader.FrameRate);
vidHeight = videoReader.Height;
vidWidth= videoReader.Width;

%% Variable initialization
labelPosition = zeros(nFrames, 2);
thresh = 30;
numFrame = 1;
xPos = zeros(nFrames,1);
yPos = zeros(nFrames,1);
flag = false;
%% Get the first frame for processing
baseR = 207;
baseG = 173;
baseB = 78;
imgWithCross = zeros(vidHeight,vidWidth,3,nFrames,'uint8');
% Extract color components
while hasFrame(videoReader)
    currrentFrame = readFrame(videoReader);
    cR = double(currrentFrame(:,:,1));
    cG = double(currrentFrame(:,:,2));
    cB =double(currrentFrame(:,:,3));
    flag = false;
    for i = 10 : vidHeight-10
        for j = 10 : vidWidth-10
            metric = sqrt((baseR-cR(i,j))^2 + (baseG-cG(i,j))^2 + (baseB-cB(i,j))^2);
            if (metric < thresh)
                xPos(numFrame) = j;
                yPos(numFrame) = i;
                flag = true;
                break;
            end 
        end
        if flag 
            break;
        end
    end
    imgWithCross(:,:,:,numFrame) = setCrosshair(currrentFrame, [yPos(numFrame) xPos(numFrame)]);
    numFrame = numFrame + 1;
end
%%
for i=1:size(imgWithCross,4)
    imshow(squeeze(imgWithCross(:,:,:,i)));
end

%%
implay(imgWithCross,30);


