function [ img ] = putCross( img,x,y )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%x=uint8(x);
%y=uint8(y);
dims=size(img);
for t=1:5%put the cross
    if rem(t,2)
        c=255;
    else
        c=0;
    end
   if (x+t<=dims(2))
       img(y,x+t)=c; 
   end
   if (x-t>=1)
       img(y,x-t)=c; 
   end
   if (y+t<=dims(1)) 
       img(y+t,x)=c; 
   end
   if (y-t>=1) 
       img(y-t,x)=c; 
   end
end
end

