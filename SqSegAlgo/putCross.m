function [ img ] = putCross( img, startPoint, crossSize)
%this function draw the cross on 3D or 2D image. 
% img - input image
% start point has to be declared as coord array: [x0 y0 z0]
startPoint=round(startPoint); % round coords if they were double 
dims=size(img);
if(numel(dims)~=numel(startPoint))
    error('Dimension mismatch! point and input image must have the same dimensions.');
end
is3D=0;
if (numel(dims)==3) % process image as 3D else 2D by default
    is3D=1; 
end
x0=startPoint(2);
y0=startPoint(1);
if(is3D)
    z0=startPoint(3);
    if(not(z0>0 && z0<=dims(3))) % if Z coord is out of range
        return;
    end
end

if(not(x0>0 && x0<=dims(1) && y0>0 && y0<=dims(2))) % if X or Y out of range
    return;
end
for t=-crossSize:crossSize%put the cross size 5 pixels 
    if rem(t,2) % choosing color black/white
        c=255; % if odd or even its distance from cross center
    else
        c=0;
    end
    if(is3D==0) 
       x=x0+t;
       if(x>0 && x<=dims(1))
           img(y0,x)=c;
       end
       y=y0+t;
       if(y>0 && y<=dims(2))
           img(y,x0)=c;
       end
    else
       x=x0+t;
       if(x>0 && x<=dims(1))
           img(y0,x,z0)=c;
       end
       y=y0+t;
       if(y>0 && y<=dims(2))
           img(y,x0,z0)=c;
       end 
       z=z0+t;
       if(z>0 && z<=dims(3))
           %img(y0,x0,z)=c;
       end 
    end
end
end

