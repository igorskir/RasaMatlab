%Формируем правильное распределение для коррекции гистограммы
p = twomodegauss(0.1, 0.05, 0.3, 0.05, 1, 0.07, 0.002);
% figure, plot(p)
% xlim([0 255])

I = imnorm(D,'norm255');
A = uint8(I);
G = histeq(A, p);
% imtool(g(:,:,20),[]);
% imtool(I(:,:,20),[]);
% imshowpair(I(:,:,20),G(:,:,20),'montage')


%СЛЕДУЮЩИЙ КОД МОЖЕТ БЫТЬ ПОЛЕЗЕН ДЛЯ ВИЗУАЛИЗАЦИИ, НО НЕ ДЛЯ ИССЛЕДОВАНИЯ
% %% STEP 3 - DISPLAYNG A 2-D CONTOUR SLICE
% %cm = brighten(jet(length(map)),-.5);
% figure
% %colormap(cm)
% contourslice(G,[],[],image_num)
% axis ij
% xlim(x)
% ylim(y)
% daspect([1,1,1])
% 
% %% STEP 4 - Displaying 3-D Contour Slices
% figure
% contourslice(G,[],[],[50],8);
% view(3);
% axis tight
