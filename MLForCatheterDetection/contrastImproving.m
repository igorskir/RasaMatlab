function [ resImage ] = contrastImproving( imgIn )
%CONTRASTIMPROVING Summary 

grayImage = imgIn;

% Compute Sobel image
[magnitudeImage, ~] = imgradient(grayImage, 'Sobel');
magnitudeImage = rangeStraching(magnitudeImage);

% Result image
resImage =  uint8(double(grayImage) - magnitudeImage);

%----block of extra code---------

% % Check that user has the Image Processing Toolbox installed.
% hasIPT = license('test', 'image_toolbox');
% if ~hasIPT
% 	% User does not have the toolbox installed.
% 	message = sprintf('Sorry, but you do not seem to have the Image Processing Toolbox.\nDo you want to try to continue anyway?');
% 	reply = questdlg(message, 'Toolbox missing', 'Yes', 'No', 'Yes');
% 	if strcmpi(reply, 'No')
% 		% User said No, so exit.
% 		return;
%     end
% end

% % Get the dimensions of the image.  
% % numberOfColorBands should be = 1.
% [rows, columns, numberOfColorBands] = size(imgIn);
% if numberOfColorBands > 1
% 	% It's not really gray scale like we expected - it's color.
% 	% Convert it to gray scale by taking only the green channel.
% 	grayImage = grayImage(:, :, 2); % Take green channel.
% end
% 

% %configurations
% format long g;
% format compact;
% fontSize = 12;
% 

% % Display the original gray scale image.
% subplot(3, 2, 1);
% imshow(grayImage, []);
% title('Original Grayscale Image', 'FontSize', fontSize);
% % Enlarge figure to full screen.
% set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
% % Give a name to the title bar.
% set(gcf,'name','Uluchshatel','numbertitle','off') 
% 
% % Compute Laplacian
% laplacianKernel = [-1,-1,-1;-1,8,-1;-1,-1,-1]/8;
% laplacianImage = imfilter(double(grayImage), laplacianKernel);
% % Display the image.
% subplot(3, 2, 2);
% imshow(laplacianImage, []);
% title('Laplacian Image', 'FontSize', fontSize);
% 
% % Compute sum of laplacian and initial image
% lapPlusInImage = double(grayImage) + laplacianImage;
% %lapPlusInImage = rangeStraching(lapPlusInImage);
% % Display the image.
% subplot(3, 2, 3);
% imshow(lapPlusInImage, []);
% title('Lap + In image ', 'FontSize', fontSize);
% 
% % Compute Sobel image
% [magnitudeImage, directionImage] = imgradient(grayImage, 'Sobel');
% % Display the gradient magnitute gray scale image.
% subplot(3, 2, 4);
% imshow(magnitudeImage, []);
% title('Sobel Image', 'FontSize', fontSize);
% 
% magnitudeImage = rangeStraching(magnitudeImage);
% % Result image
% resImage =  uint8(double(grayImage) - magnitudeImage);
% 
% resNorm = rangeStraching(resImage);
% figure; 
% imshowpair(resImage,grayImage, 'montage');
% 
% 
% % Apply an Averaging filter for Sobel image
% filteredSobelImage = conv2(magnitudeImage, ones(3)/9,'same');
% % The average Sobel
% subplot(3, 2, 5);
% imshow(filteredSobelImage, []);
% title('The average Sobel', 'FontSize', fontSize);
% 
% % Multiple the average Sobel to LaplasPlusInital
% maskImage = filteredSobelImage .* lapPlusInImage;
% 
% subplot(3, 2, 6);
% imshow(maskImage, []);
% title('maskImage', 'FontSize', fontSize);
% 
% % Result image
% maskImage = rangeStraching(maskImage);
% resImage = maskImage + double(grayImage);
% 
% figure; 
% imshowpair(resImage,grayImage, 'montage');
end

function outImg = rangeStraching(inImg)

%constant for range correction
K = 255;

minZeroImg = inImg - min(inImg(:));
maxVal = max(minZeroImg(:));

outImg = K * minZeroImg/maxVal;
end
