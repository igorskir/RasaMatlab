function [ allOk ] = CheckBorders( img, startPoint,square, direction, edgeSize, meanValue,opts)
% this function checks edges on square.
% img- input image
% start point - center of the first square
% square -array  [i j] with indexes in square Map
% direction -int 2 elements vector with values [-1/0/1, -1/0/1]
% edgeSize - currens squares edge size
center=startPoint+edgeSize*square; % coords square center on image
allOk=false;
if(direction(1)~=0 && direction(2)~=0) % Error! Propagation allowed only along X or Y axis
    return; 
end
if (direction(1)==0) % left or right square propagation
    A=center+edgeSize/2*(direction+[1, 0]); % A and B new points
    B=center+edgeSize/2*(direction+[-1, 0]);
end
if (direction(2)==0) % top or bottom propagation
    A=center+edgeSize/2*(direction+[0, 1]);
    B=center+edgeSize/2*(direction+[0, -1]);    
end

C=A-edgeSize*direction;% C and D points of already existing square
D=B-edgeSize*direction;

allOk=true;
if(ValidCoords(img,A) && ValidCoords(img,B))
    if(GetBorder(GetEdge(img,C,A),meanValue,opts)~=-1)
        allOk=false;
        return;
    end
    if(GetBorder(GetEdge(img,D,A),meanValue,opts)~=-1)
        allOk=false;
        return;
    end
    if(GetBorder(GetEdge(img,C,B),meanValue,opts)~=-1)
        allOk=false;
        return;
    end
    if(GetBorder(GetEdge(img,D,B),meanValue,opts)~=-1)
        allOk=false;
        return;
    end
    if(GetBorder(GetEdge(img,C,D),meanValue,opts)~=-1)
        allOk=false;
    end
else
    allOk=false;
end

