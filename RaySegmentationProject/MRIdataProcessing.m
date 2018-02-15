clc;
clear;
close all;
addpath(genpath(pwd))
% import square segmentation library
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

%% Setting  parameters
points=[];
dTetta=pi/2;
dFi=pi/8;

Visualize=false; % vizualize images & masks
%%
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
            startPoint(1,3) = 1; % needs for emulating 3d image 
            tic();
            %[splines, error]=SquareSegmentation(img,startPoint,options);
            error = 0; % need for sq seg. unusable for ray seg
            image3D =zeros(size(img, 1), size(img, 2), 2);
            image3D(:,:,1) = img';
            points = [];
            [points] = EmitRays(double(image3D),points,startPoint,dTetta,dFi);
            pointsForSpline = points(:,1:2); % only for xy coords
            pointsForSpline(end+1,:)= pointsForSpline(1,:);
            spline = cscvn(pointsForSpline');
            rayMask = Spline2Mask(spline,size(img));
            time=toc();
            if(error~=0)
                continue;
            end
            currAccuracy=CompareMasks(rayMask,mask);
            results(end+1,:)=[timeframe, slice, time, currAccuracy.dice];
            if(Visualize)
                imshow(img);
                figure;
                imshow(mask);
                figure;
                imshow(rayMask);
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
stdTotal = std(results(:,4));
stdTotal
failCount
imhist(results(:,4));
meanTime=mean(results(:,3));
meanTime
stdTime = std(results(:,3));
stdTime


