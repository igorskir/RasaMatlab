function [ mask ] = Splines2Mask( splines, imSize)
% returns mask builded by spline array, where 1 spline is external contour
% and others internal

externalContour=true;
for sp=1:numel(splines)
    splineContour=splines(sp); 
    t=0:0.7:splineContour.breaks(end);
    splinePoints=fnval(splineContour,t);
    x=[];
    y=[];
    for i=2:size(splinePoints,2)
        x(end+1)=splinePoints(1,i);
        y(end+1)=splinePoints(2,i);
    end
    x(end+1)=splinePoints(1,1);
    y(end+1)=splinePoints(2,1);
    if(externalContour)
        mask=poly2mask(x,y,imSize(2),imSize(1))';
        externalContour=false;
    else
        mask=mask-poly2mask(x,y,imSize(2),imSize(1))';
    end
end

end

