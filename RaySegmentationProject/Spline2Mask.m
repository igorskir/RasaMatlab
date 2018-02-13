function [ mask ] = Spline2Mask( spline, imSize)
% returns mask builded by spline
    t=0:0.7:spline.breaks(end);
    splinePoints=fnval(spline,t);
    x=[];
    y=[];
    for i=2:size(splinePoints,2)
        x(end+1)=splinePoints(1,i);
        y(end+1)=splinePoints(2,i);
    end
    x(end+1)=splinePoints(1,1);
    y(end+1)=splinePoints(2,1);
    mask=poly2mask(x,y,imSize(2),imSize(1))';
end

