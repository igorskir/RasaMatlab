function [ img ] = DrawSquare( img, center, edgeSize, color)
% draw square on image
% center - [ y x ] coords on image
edgeSize=edgeSize-1;
A=round(center+[-1 -1]*edgeSize/2);
B=round(center+[-1 1]*edgeSize/2);
C=round(center+[1 -1]*edgeSize/2);
D=round(center+[1 1]*edgeSize/2);
if(~ValidCoords(img,A) || ~ValidCoords(img,C)) % check for valid coords
    return; % not valid coords  
end

img(A(1):C(1),A(2))=color; % drawing lines of square
img(A(1),A(2):B(2))=color;
img(B(1):D(1),B(2))=color;
img(C(1),C(2):D(2))=color;

end

