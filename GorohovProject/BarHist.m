function BarHist(I)
%BARGEN generate the bar plot for image histogram 
h = imhist(I);
h10 = h(1:5:256);
horz = 1:5:256;

subplot(1,2,1);
imhist(I);
title('Subplot 1: Histogram')

subplot(1,2,2);
bar(horz, h10)
axis([0 255 0 4000])
xlabel('brighness') 
set(gca,'xtick',0:25:255)
set(gca,'ytick',0:500:4000)
title('Subplot 1: Bar Histogram')

end