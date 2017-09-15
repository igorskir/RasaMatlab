function [ ray ] = Ray(image3D,x0,y0,z0, fi, tetta, delta)
%image3D - 3D image matrix
%x0 y0 z0 - start point 
% fi and tetta - ray direction in polar coord system
% delta - step of choosing points along the ray

dx=delta*sin(tetta)*cos(fi); % projection step on each dimension 
dy=-delta*sin(tetta)*sin(fi);%
dz=delta*cos(tetta);%
dims=size(image3D); % size of input image

if(dx>0) % calculating size of output vector for each dim
   Nx=floor((dims(2)-x0)/dx);
else
    if(dx<0)
        Nx=floor(x0/abs(dx));
    else
        Nx=9999;
    end
end
% for Y axis
if(dy>0)
   Ny=floor((dims(1)-y0)/dy);
else
    if(dy<0)
        Ny=floor(y0/abs(dy));
    else
        Ny=99999;
    end
end
% for Z axis
if(dz>0)
   Nz=floor((dims(3)-z0)/dz);
else
    if(dz<0)
        Nz=floor(z0/abs(dz));
    else
        Nz=99999;
    end
end

N=min([Nx Ny Nz]);% choosing the minimal length from all dims
ray=zeros(N,1); % creating output ray

for i=0:N-1
    x=round(x0+dx*i);
    y=round(y0+dy*i);
    z=round(z0+dz*i);
    ray(i+1)=image3D(y,x,z);
end

