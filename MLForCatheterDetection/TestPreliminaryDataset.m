ax = 'long1';
minSlice = 55;
maxSlice = 115;
for nTimeframe = 1:17
    filename = 'LV Catheter 07.nrrd'; % delete after getting  features
    [X, meta] = nrrdread(filename);
    sz = sscanf(meta.sizes, '%d');
    nDims = sscanf(meta.dimension, '%d');
    I = squeeze(X(:,:,:,nTimeframe));
    
    for nSlice = minSlice:maxSlice
        img = GetImage(I, nSlice, ax);
        [level,EM] = graythresh(img);
        BW = imbinarize(img, level);
        hFig = figure;      
        imshowpair(img, BW, 'montage');
        str0 = sprintf("%d timeframe: %d:%d", nTimeframe, minSlice, maxSlice);
        str1 = sprintf("%d slice", nSlice);
        AddTitle({str0, str1});
        set(gcf, 'Position', [scrSz(3), 0, scrSz(3), scrSz(4)],...
        'Color', 'w', 'name', str1, 'numbertitle', 'off'); % FOR THE SECOND DISPLAY ONLY
        close(hFig);
    end
end


