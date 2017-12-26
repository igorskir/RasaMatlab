function [ squaresNew, mapNew, currId, procListCount ] = CrushSquares(map,mapCenter,currId)
% this function crush each square to 4 squares
[rows, cols]=size(map);
squaresNew=zeros(rows*cols,2);
mapNew=zeros(rows,cols);
Ic=mapCenter(1);
Jc=mapCenter(2); % extracting center coords
procListCount=0;
for i=1:rows
    for j=1:cols
        if (map(i,j)>0) % if cell is a square crushing it...
            Inew=Ic+2*(i-Ic); % stretching  coords
            Jnew=Jc+2*(j-Jc);
            mapNew(Inew,Jnew)=map(i,j);
            
            % marking uo 3 additional squares
            mapNew(Inew,Jnew+1)=currId;
            mapNew(Inew+1,Jnew)=currId+1;
            mapNew(Inew+1,Jnew+1)=currId+2;
            currId=currId+3;
            
            %borders for current square
            topBorder=map(i-1,j)>0; 
            bottomBorder=map(i+1,j)>0;
            leftBorder=map(i,j-1)>0;
            rightBorder=map(i,j+1)>0; 
            % adding outside squares to processing list
            if(~topBorder || ~leftBorder)
                squaresNew(procListCount+1,:)=[Inew, Jnew]-mapCenter;
                procListCount=procListCount+1; 
            end
            if(~topBorder || ~rightBorder)
                squaresNew(procListCount+1,:)=[Inew, Jnew+1]-mapCenter;
                procListCount=procListCount+1; 
            end
            if(~bottomBorder || ~leftBorder)
                squaresNew(procListCount+1,:)=[Inew+1, Jnew]-mapCenter;
                procListCount=procListCount+1; 
            end
            if(~bottomBorder || ~rightBorder)
                squaresNew(procListCount+1,:)=[Inew+1, Jnew+1]-mapCenter;
                procListCount=procListCount+1; 
            end
        end
    end
end
    
end

