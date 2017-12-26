function [ valid ] = ValidCoords( img, point )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
valid=1;
point=round(point);
[N M]=size(img);
if(point(1)<=0 || point(2)<=0 || point(1)>N || point(2)>M)
    valid=0;
end

end

