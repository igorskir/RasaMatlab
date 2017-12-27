%% Test one-feature accuracy
model = zeros(1,20);
filename = 'Feature engineering.xlsm';
sheet = 'Catheter analysis';
row = 115;
column = 'D';

for idx = 1:20
    model(1,idx) = 1;
    NetTrain;
    model = zeros(1,20);
    xlRange = strcat(num2str(row), char(column+idx-1));
    xlswrite(filename, cathClassRate, sheet, xlRange)
    fprintf('%d iteration \n', idx);
end
