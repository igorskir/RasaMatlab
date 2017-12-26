function [ pp ] = GetSplineContours( img, startPoint, map, mapCenter, edgeSize, meanValue,opts)
% returns piecewise polynoms which describes contours 
% first -external
% others - internal
pp=repmat(cscvn([1 1; 1 2]),1,5); % creating array of
currPolyIndex=1;
isExternal=true;
for i=1:size(map,1)
    for j=1:size(map,2)
        if(map(i,j)==-1)
            % check if it's not one cell area
            dir=[1 0];
            oneCellArea=true;
            for direction=1:4
                if(map(i+dir(1),j+dir(2))<=0)
                    oneCellArea=false;
                end
                dir=TurnLeft(dir);
            end
            if(oneCellArea)
                continue;
            end
            
            if(isExternal) % first contour is always external
                isExternal=false;
                dir=[0 1]; % direcction corresponds clockwise order
            else % for internal areas counter-clockwise
                dir=[0 -1];
                if(map(i-dir(1),j-dir(2))>0)
                    dir=TurnRight(dir);
                end
            end
            currPoint=1; % pointer  to the first free row in point list
            points=zeros(40000,2); % !!! CRUTCH!!!!
            startSq=[i j];
            currSq=startSq;
            endCondition=true;
            while(endCondition)
                sqCenter=startPoint+edgeSize*(currSq-mapCenter);
                % looking for a border...
                A=sqCenter-edgeSize/2*dir+edgeSize/2*TurnRight(dir);
                B=sqCenter-edgeSize/2*dir+edgeSize/2*TurnLeft(dir);
                index=GetBorder(GetEdge(img,A,B),meanValue,opts);
                if(index~=-1)
                    points(currPoint,:)=A+(index-1)*TurnLeft(dir); 
                    currPoint=currPoint+1; 
                end
                if(index==-1 && ~ValidCoords(img,B))
                  points(currPoint,:)=CorrectNotValidPoint(img,B); 
                  currPoint=currPoint+1;   
                end
                % turn letft while cant go forward
                while(map(currSq(1)+dir(1),currSq(2)+dir(2))>0)
                    dir=TurnLeft(dir);
                end
                nextSq=currSq+dir+TurnRight(dir);
                if(map(nextSq(1),nextSq(2))<0)
                    A=sqCenter+edgeSize/2*dir+edgeSize/2*TurnRight(dir);
                    B=sqCenter+edgeSize/2*dir+edgeSize/2*TurnLeft(dir);
                    index=GetBorder(GetEdge(img,A,B),meanValue,opts);
                    if(index~=-1)
                        points(currPoint,:)=A+(index-1)*TurnLeft(dir); 
                        currPoint=currPoint+1;
                    end
                    if(index==-1 && ~ValidCoords(img,B))
                      points(currPoint,:)=CorrectNotValidPoint(img,B); 
                      currPoint=currPoint+1;   
                    end
                    map(currSq(1),currSq(2))=-2;
                    currSq=nextSq;
                    dir=TurnRight(dir);
                else
                    map(currSq(1),currSq(2))=-2;
                    currSq=currSq+dir;
                end
                endCondition=~(currSq(1)==startSq(1) && currSq(2)==startSq(2));
            end
            if(currPoint~=1)
                points(currPoint,:)=points(1,:); 
                points=points(1:currPoint,:); % removing empty rows
                pp(currPolyIndex)=cscvn(points');
                currPolyIndex=currPolyIndex+1;
            end
        end
    end
end
pp=pp(1:currPolyIndex-1); % removing empty elemets
end

