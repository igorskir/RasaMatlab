clc;
clear;
% import square segmentation library
addpath('SqSegAlgo');
% loading sol data
fileNum=11;
images=load(['YORK\Raw_data\sol_yxzt_pat' num2str(fileNum) '.mat']);
images=images.sol_yxzt;
contours=load(['YORK\Segmented_data\manual_seg_32points_pat' num2str(fileNum) '.mat']);
contours=contours.manual_seg_32points;
startPoints=load(['startPoints' num2str(fileNum) '.mat']);
startPoints=startPoints.startPoints;

[height, width, ~, ~]=size(images);
masks=GetMasks(contours,[height, width]);

% SqSeg options
options.threshold=21.2; % 55.5 optimal for twomedegauss
options.level=6.2; % 12.1 optimal
options.edgeSize=10; % 10-5 optimal
options.minEdgeSize=5;
Visualize=false; % vizualize images & masks
results=[];
for slice=1:size(masks,1)
    for timeframe=1:size(masks,2)
        img=uint8(squeeze(images(:,:,slice,timeframe)));
        mask=masks{slice,timeframe};
        %imshow(img);
        if (mask~=-99999)
            % starting test
            p = twomodegauss(0.1, 0.05, 0.3, 0.05, 1, 0.07, 0.002);
            %imhist(img);
            %figure;
            img=histeq(img,255);
            %imhist(img);
            %figure;
            %imshow(img);
            %return;
            startPoint= startPoints{slice,timeframe};
            tic();
            [splines, error]=SquareSegmentation(img,startPoint,options);
            time=toc();
            if(error~=0)
                continue;
            end
            sqMask=Splines2Mask(splines,size(img));
            currAccuracy=CompareMasks(sqMask,mask);
            results(end+1,:)=[timeframe, slice, time, currAccuracy.dice];
            if(Visualize)
                imshow(img);
                figure;
                %imshow(mask);
                %figure;
                imshow(sqMask);
                %ginput(1);
                close all;
            end
        end
    end
end
total=0;
totalCount=0;
failCount=0;
for i=1:size(results,1)
    if(results(i,4)>0.5)
        total=total+results(i,4);
        totalCount=totalCount+1;
    else
        failCount=failCount+1;
    end
end
total=total/totalCount;
total
failCount
imhist(results(:,4));
meanTime=mean(results(:,3));
meanTime


