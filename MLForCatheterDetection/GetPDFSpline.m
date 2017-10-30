function [ pp ] = GetPDFSpline(data,binCount)
maxEl=max(data);
minEl=min(data);
%binCount=10;
step=(maxEl-minEl)/binCount;
edges=minEl:step:maxEl;
bins=zeros(binCount,1);

for i=1:numel(data)
    for j=1:binCount
        currEl=data(i);
        if(currEl>=edges(j) && currEl<=edges(j+1))
            bins(j)=bins(j)+1;
            break;
        end
    end
end
bins=bins/(numel(data)*step);
middles=zeros(binCount,1);
for i=1:binCount
    middles(i)=(edges(i)+edges(i+1))/2;
end
figure;
plot(middles,bins);
hold on;
x=minEl:step*0.01:maxEl;
y=spline(middles,bins,x);
sum(y)*step*0.01
plot(x,y,'red');
pp=cscvn([middles; bins]);
end
