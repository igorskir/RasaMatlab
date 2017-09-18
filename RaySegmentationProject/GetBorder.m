function [ index ] = GetBorder( ray, delta, width,level ) 
%looking for the the border on a ray
%width - width of approximation window 
ray=smooth(ray,4);
%plot(ray);
i=width+1;
x=0:delta:(2*width)*delta;
k=0;
b=0;
flag=0;
while((i<=(numel(ray)-width)) && flag==0)
    y=ray(i-width:i+width);
    %coeffs=[1 1];
    coeffs=polyfit(x',y,1); %coeffs of linear approximation
    i=i+1;
    k=abs(coeffs(1));
    %k
    if(k>=level) 
        flag=1; 
    end %threshold level;
end
index=i;
if(flag==0) % if cycle ends without results 
   index=-1;
end
end

