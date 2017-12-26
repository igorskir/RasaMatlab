function [ outputSplines, errorCode] = SquareSegmentation(img,startPoint,opts)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
errorCode=0;
edgeSize=opts.edgeSize;
minEdgeSize=opts.minEdgeSize;
assert(floor(log2(edgeSize/minEdgeSize))-log2(edgeSize/minEdgeSize)==0,'Incompatable edge sizes!')
%putting first square
directions=[[0 1]; [-1 0]; [0 -1]; [1 0]];% setting up parameters
allDiagonalsCorrect=true;
img=double(img);
sum=img(round(startPoint(1)),round(startPoint(2))); % sum of brightnesses for mean value

brCount=1; %  count of points for mean brightness
for dx=-1:2:1
    for dy=-1:2:1 % checking diagonals
            endPoint=CorrectNotValidPoint(img,startPoint+[dx dy]*edgeSize/2);
            sum=sum+img(round(endPoint(1)),round(endPoint(2))); % mean brightness: adding new point
            brCount=brCount+1;
            if (GetBorder(GetEdge(img,startPoint,endPoint),sum/brCount,opts)~=-1)
                allDiagonalsCorrect=false;
            end
    end
end
%assert(allDiagonalsCorrect, 'Impossible to put first square! Try to change edge size or start point!');
if(~allDiagonalsCorrect) % error! stop executing...
    errorCode=1;
    outputSplines=0;
    return;
end
[height, width]=size(img); % get image size
map=zeros(floor(height/minEdgeSize)+20,floor(width/minEdgeSize)+20); 
%map of squares: each cell contains square Id or null
mapCenter=[round(startPoint(1)/minEdgeSize)+10 round(startPoint(2)/minEdgeSize)+10];
squares=zeros(size(map,1)*size(map,2),2);
squares(1,:)=[0 0];% matrix: each row is X and Y indexes of square in map
procListCount=1; % Length of processing list squares

currId=1; % current identifier 
map(mapCenter(1),mapCenter(2))=currId;% mark up first square on the map
currId=currId+1; % choosing next free Id for square
while (edgeSize>=minEdgeSize)
    currSquare=1; % starts from first square in list
    borderSquaresCount=0;
    while(currSquare<=procListCount) 
        square=squeeze(squares(currSquare,:));
        for direction=1:4 % trying to create a neighbour sq for each side
           neighb=mapCenter+square+squeeze(directions(direction,:)); % getting neighbour sq
           if(~ValidCoords(map,neighb)) % if neighbour out of bounds
               continue; % 
           end
            if(map(neighb(1),neighb(2))==0) % side is free
                  if(CheckBorders(img,startPoint,neighb-mapCenter,squeeze(directions(direction,:)),edgeSize,sum/brCount,opts))
                       squares(procListCount+1,:)=neighb-mapCenter;
                       procListCount=procListCount+1;
                       map(neighb(1),neighb(2))=currId;
                       currId=currId+1;
                       % mean brightness of region correctin
                       point=startPoint+edgeSize*(neighb-mapCenter);
                       sum=sum+img(round(point(1)),round(point(2))); % mean brightness: adding new point
                       brCount=brCount+1;
                  else
                      map(neighb(1),neighb(2))=-1;
                      borderSquaresCount=borderSquaresCount+1;
                  end
            end
        end
        currSquare=currSquare+1;
    end
    % ------------------------------------------ DEBUG DRAWING
%     figure;
%     imgForDraw=uint8(img);
%     for i=1:size(map,1)
%         for j=1:size(map,2)
%             if(map(i,j)~=0)
%                     if(map(i,j)<0)
%                         color=150;
%                     else
%                         color=0;
%                     end
%                     imgForDraw=DrawSquare(imgForDraw,startPoint+edgeSize*([i,j]-mapCenter),edgeSize,color); % let's draw it!
%             end
%         end
%     end
%     imshow(imgForDraw);
    % ---------------------------------------- END DEBUG DRAWING
    if(edgeSize<=minEdgeSize)
        break;
    end
    % changing edge size and crushing squares
    [squares, map, currId, procListCount]=CrushSquares(map,mapCenter,currId);
    edgeSize=edgeSize/2;
    startPoint=startPoint-[edgeSize/2, edgeSize/2];
end

outputSplines=GetSplineContours( img, startPoint, map, mapCenter, edgeSize, sum/brCount,opts);

