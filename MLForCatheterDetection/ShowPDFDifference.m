function [ ] = ShowPDFDifference( dataTissue, dataCath, binCountTissue, binCountCath)

[xTissue, probTissue]=GetPDFSpline(dataTissue,binCountTissue);
[xCath, probCath]=GetPDFSpline(dataCath,binCountCath);


step=(xCath(end)-xCath(1))/(binCountCath*20);
sum=0;
tissueMin=min(xTissue);
tissueMax=max(xTissue);
for i=xCath(1):step:xCath(end)
    splineCath=spline(xCath, probCath,i);
    if (i>=tissueMin && i<=tissueMax)
        splineTissue=spline(xTissue, probTissue,i);
        
        if (splineTissue<=splineCath)
            sum=sum+(splineCath-splineTissue)*step;
        end
    else
        sum=sum+splineCath*step;
    end
end
plot(xTissue, probTissue);
hold on;
plot(xCath, probCath);
grid on;
sum

end

