function [middleVals, bins] = GetPDFSpline(data, binCount, isDiscrete)
maxElement = max(data);
minElement = min(data);
middleVals = 0;

if ~isDiscrete 
    data=sort(data);
    elementsPerBin=floor(numel(data)/binCount)+1;
    %binCount=floor(numel(data)/elementsPerBin)+1;
    %edges=zeros(binCount+1,1);
    %bins=zeros(binCount,1);
    edges(1)=data(1);
    index=elementsPerBin;
    i=1;
    while(index<numel(data))
            bins(i)=elementsPerBin;
            edges(i+1)=(data(index)+data(index+1))/2;
            index=index+elementsPerBin;
            i=i+1;
    end
    %edges(end)=data(end);
    index=index-elementsPerBin;
    bins(i)=(numel(data)-index);
    edges(i+1)=data(end);
    binCount=numel(bins);
    % Merging bad bins
    allMerged=false;
    while (~allMerged)
        allMerged=true;
        for i=1:binCount
            if(edges(i)==edges(i+1))
                allMerged=false;
                if(i==1) % right merge
                    mergeType='R';
                else
                    if(i==binCount) % left merge
                        mergeType='L';
                    else
                        d1=edges(i)-edges(i-1);
                        d2=edges(i+1)-edges(i);
                        if(d1<d2)
                            mergeType='L';
                        else
                            mergeType='R';
                        end
                    end
                end
                if(mergeType=='L')
                    bins(i-1)=bins(i-1)+bins(i);
                    bins(i)=[];
                    edges(i)=[];
                else
                    bins(i+1)=bins(i+1)+bins(i);
                    bins(i)=[];
                    edges(i+1)=[];
                end
                binCount=binCount-1;                
                break;
            end
        end      
    end
    middleVals=zeros(binCount,1);
    for i=1:binCount
        middleVals(i)=(edges(i)+edges(i+1))/2;
        floatStep=edges(i+1)-edges(i);
        bins(i)=bins(i)/(numel(data)*floatStep); % comment for the StemPlot
    end
%     bins(end+1) = bins(end); % uncomment for the StemPlot
%     StemPlot(edges, bins);   % uncomment for the StemPlot

%     step = (maxElement - minElement)/binCount;
%     edges = minElement:step:maxElement;
%     bins = zeros(binCount, 1);
%     for i = 1:numel(data)
%         for j = 1:binCount
%             currElement = data(i);
%             if(currElement >= edges(j) && currElement <= edges(j+1))
%                 bins(j) = bins(j)+1;
%                 break;
%             end
%         end
%     end
%     bins = bins/(numel(data)*step);
%     middleVals = zeros(binCount, 1);
%     for i = 1:binCount
%         middleVals(i) =(edges(i)+edges(i+1))/2;
%     end
else
   binCount = maxElement - minElement + 1;
   edges = binCount;
   bins = zeros(binCount, 1);
       for i = 1:numel(data)      
            currElement = data(i);
            bins(currElement - minElement + 1) =  bins(currElement - minElement + 1) + 1;
       end
   bins = bins/numel(data);
   middleVals=minElement:1:maxElement;
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
