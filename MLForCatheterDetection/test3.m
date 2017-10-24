%% Extract relevant input feature array (manual SFS)
filename = 'D:\Clouds\Google Drive\RASA\Matlab\12. Catheter ML\Presentation\Feature engineering.xlsm';
sheet = 5;
rng = 'all'; % all - 18, big - 11, mid - 9, small - 5
switch rng
    case 'all'
        xlRange = 'D110:U110';
    case 'big'
        xlRange = 'D111:U111';  
    case 'mid'
        xlRange = 'D112:U112';
    case 'small'
        xlRange = 'D113:U113';
end

model = xlsread(filename,sheet,xlRange);
[numRows, numCols] = size(netTrainInputs);
netSelectInputs = zeros(numRows, 1);
for i = 1:numCols
        if model(1, i) == 1
            netSelectInputs = [netSelectInputs, x(1:end, i)];
        end
end
netSelectInputs(:,1) = [];
vars.featureSelection = {'filename', 'sheet', 'rng', 'model', 'numRows',...
                         'numCols','i'};
clear(vars.featureSelection{:});
