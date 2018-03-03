function output = GetSeparatedRegions(BW, method, ax)
% se = strel('disk', 2); %
% se = strel('ball', 2, 0, 0);
D = -bwdist(~BW, method);
%imshow(D, 'InitialMagnification', 'fit');
    if strcmp(ax, 'short')   
        D(~BW) = Inf;
        L = watershed(D);
        %imshow(label2rgb(L, colormap(brewermap([],'Set3'))), 'InitialMagnification', 'fit')
        L(~BW) = 0;
        %imshow(label2rgb(L, colormap(brewermap([],'Set3'))), 'InitialMagnification', 'fit')
        output = logical(L);
    elseif strcmp(ax, 'long') || strcmp(ax, 'long1') || strcmp(ax, 'long2')
        mask = imextendedmin(D, 1);
        % imshowpair(BW,mask,'blend')
        D = imimposemin(D, mask);
        L = watershed(D);
        output = BW;
        output(L == 0) = 0;
    end
% output = imopen(output, se);
end

