function [ edge ] = GetEdge( img, startPoint, endPoint)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
distance=sqrt(sum((startPoint-endPoint).^2));
len=round(distance);
edge=uint8(zeros(len+1+4,1));
direction=(endPoint-startPoint)/len;
point=startPoint-2*direction;
for i=-2:len+2
    pCorrect=CorrectNotValidPoint(img,point);
    edge(i+3)=img(round(pCorrect(1)),round(pCorrect(2)));
    point=point+direction;
end
% plot(edge); % DEBUG PLOT!
% figure;
% im2=putCross(img,startPoint,2);
% im2=putCross(im2,endPoint,1);
% imshow(im2);
end

