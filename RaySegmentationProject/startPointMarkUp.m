% this script marks up start points for SqSegAlgo and save them as cell
% array
clc;
clear;
fileNum=11; %loading images and contours
images=load(['YORK\Raw_data\sol_yxzt_pat' num2str(fileNum) '.mat']);
images=images.sol_yxzt;
contours=load(['YORK\Segmented_data\manual_seg_32points_pat' num2str(fileNum) '.mat']);
contours=contours.manual_seg_32points;
% cell array for coords
startPoints=cell(size(images,3),size(images,4));
half=34; % half o the set. el #33 separate two contours: external and internal
for slice=1:size(images,3)
    for timeframe=1:size(images,4)
        img=uint8(squeeze(images(:,:,slice,timeframe)));
        if (contours{slice,timeframe}~=-99999)
            img=histeq(img);
            figure;
            imshow(img);
            contour=contours{slice,timeframe}(1:half-2,:);
            for t=2:size(contour,1) % ploting contour over current image
                l=line([contour(t-1,1) contour(t,1)], [contour(t-1,2) contour(t,2)]);
                set(l, {'color','LineWidth'}, {'red',1});
            end
            l=line([contour(end,1) contour(1,1)], [contour(end,2) contour(1,2)]);
            set(l, {'color','LineWidth'}, {'red',1});
            coords=ginput(1); % get mouse click coords
            startPoints{slice,timeframe}=squeeze( [coords(1,2), coords(1,1)]);  
            close all;
        else
            startPoints{slice,timeframe}=-99999;
        end
    end
end
save(['startPoints' num2str(fileNum) '.mat'],'startPoints');