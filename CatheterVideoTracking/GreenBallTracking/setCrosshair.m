function img = setCrosshair( img, position )
sz = 7;
for i = -sz:sz
    for j = -sz:sz
        curPos = position + [i j];
        img(curPos(1), curPos(2),  1) = uint8(255);
        img(curPos(1), curPos(2),  2) = uint8(0);
        img(curPos(1), curPos(2),  3) = uint8(0);
    end
end
end

