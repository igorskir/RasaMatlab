function img = GetImage(I, nSlice, ax)

    switch ax
        case 'short'
            img = I(:,:,nSlice);
        case 'long1'
            tempImg = squeeze(I(nSlice,:,:));
            img = rot90(tempImg);
        case 'long2'
            tempImg = squeeze(I(:,nSlice,:));
            img = rot90(tempImg);
        otherwise
            disp('Wrong axis value')
    end

end

