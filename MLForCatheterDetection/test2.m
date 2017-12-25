%% Reading the data
tic;
filename = 'LV Catheter 07.nrrd'; % delete after getting  features
[X, meta] = nrrdread(filename);
Y = double(X);
sz = sscanf(meta.sizes, '%d');
nDims = sscanf(meta.dimension, '%d');
toc;
%%
I = load('I.mat');
I = I.I;
%%
numCases = zeros(1,1);
diceCases = zeros(1,1);
j = 0;
for nSlice = 50:110 
img = I(:,:,nSlice);
sz = size(img);
level = graythresh(img);


BW = imbinarize(img, level);
[BWmser, count] = DetectMSERblobs(BW);
numCases(nSlice-49,1) = count; 

if count ~= 1
    continue;
end
BWref = RealCatheterDrawing(BWmser);
idxDice = dice(BWref, BWmser);
if idxDice ~= 0
    j = j + 1; 
end
diceCases(j,1) = idxDice; 
end
