function [ img ] = putCross( img,x,y )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

x=round(x);
y=round(y);
if(x<=0 || y<=0 || x>size(img,2) || y>size(img,1))
    return;
end;
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

