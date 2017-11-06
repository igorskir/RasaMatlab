for i = 30:50 % 35
    img = dicomread('HAV6GEG2.dcm');
    gray = rgb2gray(img(:,:,:,numSlice));
    [level,EM] = graythresh(gray);
    BW = imbinarize(gray, level);
    hFig = figure;
    imshowpair(gray, BW, 'montage');
    close(hFig)
end