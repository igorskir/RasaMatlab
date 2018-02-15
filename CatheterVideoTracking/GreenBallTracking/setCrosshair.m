function img = setCrosshair( img, position )
position = double(position);
sz = 7;
for i = -sz:sz
    for j = -sz:sz
        curPos = position + [i j];
        img(curPos(1), curPos(2),  1) = 255;
        img(curPos(1), curPos(2),  2) = 0;
        img(curPos(1), curPos(2),  3) = 0;
    end
end
end

