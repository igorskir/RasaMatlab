function [middleVals, bins] = GetPDFSpline(data, binCount)
maxElement = max(data);
minElement = min(data);
step = (maxElement - minElement)/binCount;
edges = minElement:step:maxElement;
bins = zeros(binCount, 1);
for i = 1:numel(data)
    for j = 1:binCount
        currElement = data(i);
        if(currElement >= edges(j) && currElement <= edges(j+1))
            bins(j) = bins(j)+1;
            break;
        end
    end
end
bins = bins/(numel(data)*step);
middleVals = zeros(binCount, 1);
for i = 1:binCount
    middleVals(i) =(edges(i)+edges(i+1))/2;
end
% figure;
% plot(middles,bins,'*');
% hold on;
% x=minEl:step*0.01:maxEl;
% %y=spline(middles,bins,x);
% %sum(y)*step*0.01
% %plot(x,y,'red');
% grid on
% pp=cscvn([]);
% y=ppval(pp,x);
% y
% hold on;
% plot(x,y);
end
