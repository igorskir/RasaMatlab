function [ point ] = CorrectNotValidPoint( img, point)
[N, M]=size(img);
if(point(1)<1)
    point(1)=1;
end
if(point(2)<1)
    point(2)=1;
end
if(point(1)>N)
    point(1)=N;
end
if(point(2)>M)
    point(2)=M;
end
end