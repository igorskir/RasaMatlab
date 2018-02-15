function [ points ] = EmitRays( image3D, points, startPoint, dTetta, dFi)
% this function emulate ray emission from startPoint [x0 y0 z0]
% in spherical coords using steps dTetta and dFi for each coord
%  exceptDirection - direction [tetta fi] which should not be considered
%   points - in this list will be added detected points

%  optimization mark: for tetta=0 and pi its enough to get Ray once without
%  looping on fi. UPDATED
x0=startPoint(1); % extracting values...
y0=startPoint(2);
z0=startPoint(3);
%params for border detection
level=4.5; % condition for linear approximation coeffs
width=2; % width of a window for calculating approximation coeffs
delta=1;
maxDist=12.3;

for tetta=pi/2:dTetta:pi/2
    for fi=0:dFi:(2*pi-dFi)  
        % here we need to check current [tetta fi] if its not excepted direction
        % not released yet
        ray=Ray(image3D, startPoint,fi,tetta,delta); % creating ray for current direction
        borderIndex=GetBorder(ray,delta,width,level,maxDist); % looking for a border
        if(borderIndex==-1)
               continue;
        end
        dx=delta*sin(tetta)*cos(fi); % projection step on each dimension 
        dy=-delta*sin(tetta)*sin(fi);%
        dz=delta*cos(tetta);%
        x=round(x0+dx*borderIndex);
        y=round(y0+dy*borderIndex);
        z=round(z0+dz*borderIndex);
        points(size(points,1)+1,:)=[x y z image3D(x,y,z)];
        if(tetta==0 || tetta==pi) % for this values enough first iteration
            break;
        end
    end
end
% avg=mean(squeeze(points(:,4))); % calculating average intensity
% dev=std(squeeze(points(:,4))); % and standard deviation for border points
% i=1;
% while (i<=size(points,1)) % drop out points with intensity out of range [avg-3*dev .. avg+3*dev]
%     if(abs(points(i,4)-avg)>=dev) 
%         points(i,:)=[]; 
%     else
%         i=i+1;
%     end
% end
pointCount=size(points,1);
distances=zeros(pointCount); %matrix of distances between points. distance(i,j)-distance between I-th and J-th
for i=1:pointCount
    for j=1:i
        if(i==j)
            distances(i,j)=9999;
            continue;
        end
        distances(i,j)=0;
        for d=1:3 % for each dimension XYZ
            distances(i,j)=distances(i,j)+((points(i,d)-points(j,d))^2);
        end
        distances(i,j)=sqrt(distances(i,j));
        distances(j,i)=distances(i,j); % symmetry mapping
    end
end
% [maxDist, maxIndex]=max(min(distances)); % detecting the furthest point
%maxDist 
% BELOW: calculations of new start point coords and calling  THIS func
% recursively 
% nextCenter=[(x+points(maxIndex,1))/2, (y+points(maxIndex,2))/2, (z+points(maxIndex,3))/2];
% if(maxDist>35)
%     EmitRays(image3D,points,nextCenter,dTetta,dFi,[]);
% end
end

