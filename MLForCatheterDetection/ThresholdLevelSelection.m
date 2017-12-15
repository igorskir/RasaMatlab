%% Reading the data
isVisual = 0;
nTimeframe = 9; %9
scrSz = get(0, 'Screensize');
tic;
filename = 'LV Catheter 07.nrrd'; % delete after getting  features
[X, meta] = nrrdread(filename);
Y = double(X);
sz = sscanf(meta.sizes, '%d');
nDims = sscanf(meta.dimension, '%d');
toc;

%% Estimation
diceDataset = [];
for nTimeframe = 1:17
    % I = load('I.mat');
    % I = I.I;
    I = squeeze(X(:,:,:,nTimeframe));
    %%
    numCases = zeros(1,1);
    diceTimeframe = zeros(1,1);
    j = 0;
    for nSlice = 69:105 
        img = I(:,:,nSlice);
        sz = size(img);
        level = GetThresholdingLevel(img, 'OtsuHistogram');
        BW = imbinarize(img, level);
        [BWmser, count] = DetectMSERblobs(BW);
        numCases(nSlice-58,1) = count; 
        if count ~= 1
            continue
        end
        BWref = RealCatheterDrawing(BWmser);
        idxDice = dice(BWref, BWmser);
        if isVisual == 1
            C = imfuse(BWref,BWmser, 'ColorChannels', [1 2 0]);
            % imshow(C, 'InitialMagnification', 'Fit')
            imshowpair(img, C, 'montage');
            set(gcf, 'Position', [scrSz(3), 0, scrSz(3), scrSz(4)],...
                'Color', 'w', 'numbertitle', 'off');
            str = sprintf("Dice index: %.2f", idxDice);
            AddTitle(str);
        end

        if idxDice ~= 0
            j = j + 1; 
        end
        diceTimeframe(j,1) = idxDice; 
    end
    diceDataset = vertcat(diceDataset, diceTimeframe);
end