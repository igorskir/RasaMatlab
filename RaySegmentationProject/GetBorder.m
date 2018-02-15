function [ index ] = GetBorder( ray, delta, width,level, maxDist) 
%looking for the the border on a ray
%width - width of approximation window 
raw=ray;
%ray=smooth(ray,10);
width = 2;
%plot(ray);
i=width+1;
x=0:delta:(2*width)*delta;
b=0;
s=sum(raw(1:i-1));
while((i<=(numel(ray)-width)))
    y=ray(i-width:i+width);
    k=LSM(y');    
    if(k>=level || abs(s/(i-1)-raw(i))>maxDist) 
        break; 
    end %threshold level;
    s=s+raw(i);
    i=i+1;
end
index=i;
if(flag==0) % if cycle ends without results 
   index=-1;
end
end

