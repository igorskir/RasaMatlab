function [ masks ] = GetMasks(contours, imSize)

half=34;
masks=cell(size(contours));
for timeframe=1:size(contours,2)
    for zVal=1:size(contours,1)
        if(contours{zVal,timeframe}~=-99999) % processing...
           nodes=contours{zVal,timeframe};
           %xNodes=nodes(half:end,2);
           %yNodes=nodes(half:end,1);
           xNodes=nodes(1:half-2,2);
           yNodes=nodes(1:half-2,1);
           xNodes(end+1)=xNodes(1); 
           yNodes(end+1)=yNodes(1);
           masks{zVal,timeframe}=poly2mask(xNodes,yNodes,imSize(2),imSize(1))';
        else
           masks{zVal,timeframe}=-99999; 
        end
    end

end

