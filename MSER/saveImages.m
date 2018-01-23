nSlice = 53:108;
for i = nSlice
   tempImg = BWr(:,:,i);
   imgName = [num2str(i),'.png'];
   imwrite(tempImg, imgName);
end

