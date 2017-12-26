function [ index ] = GetBorder( ray, meanValue, opts) 
% 1D region growing 
%looking for the the border on a ray/edge
ray=double(ray);

maxDist=opts.threshold;
kLimit=opts.level;

index=-1;
i=3;

while(i<=numel(ray)-2) 
    [k]=LSM(ray(i-2:i+2)');
    if((abs(ray(i)-meanValue)>maxDist) || (abs(k)>kLimit))
        index=i-2;
        break;
    end
    i=i+1;
end
end

