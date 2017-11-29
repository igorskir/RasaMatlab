function [ imgOut ] = contrastImproving( imgIn )
%CONTRASTIMPROVING Summary of this function goes here
%   Detailed explanation goes here

%configurations
format long g;
format compact;
fontSize = 12;

% Check that user has the Image Processing Toolbox installed.
hasIPT = license('test', 'image_toolbox');
if ~hasIPT
	% User does not have the toolbox installed.
	message = sprintf('Sorry, but you do not seem to have the Image Processing Toolbox.\nDo you want to try to continue anyway?');
	reply = questdlg(message, 'Toolbox missing', 'Yes', 'No', 'Yes');
	if strcmpi(reply, 'No')
		% User said No, so exit.
		return;
    end
end

grayImage = imgIn;
% Get the dimensions of the image.  
% numberOfColorBands should be = 1.
[rows, columns, numberOfColorBands] = size(imgIn);
if numberOfColorBands > 1
	% It's not really gray scale like we expected - it's color.
	% Convert it to gray scale by taking only the green channel.
	grayImage = grayImage(:, :, 2); % Take green channel.
end

% Display the original gray scale image.
subplot(2, 2, 1);
imshow(grayImage, []);
title('Original Grayscale Image', 'FontSize', fontSize);
% Enlarge figure to full screen.
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
% Give a name to the title bar.
set(gcf,'name','Demo by ImageAnalyst','numbertitle','off') 

% Compute Laplacian
laplacianKernel = [-1,-1,-1;-1,8,-1;-1,-1,-1]/8;
laplacianImage = imfilter(double(grayImage), laplacianKernel);
% Display the image.
subplot(2, 2, 2);
imshow(laplacianImage, []);
title('Laplacian Image', 'FontSize', fontSize);

sharpenedImage = double(grayImage) + laplacianImage;
% Display the image.
subplot(2, 2, 3);
imshow(sharpenedImage, []);
title('Sharpened Image', 'FontSize', fontSize);

end

function outI = helper(img)

minI = min(img(:));
K = 255;
fm = img - minI;
maxI = max(fm(:));

fs =fm/maxI;
outI = K * fs;

end
