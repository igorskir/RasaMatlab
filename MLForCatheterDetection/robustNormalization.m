function output = robustNormalization(input, normType, isNormFirstCol)

if isstruct(input)
    [inputArray, fields] = struct2mat(input);
    fields = fields';
    numRows = numel(input);
    numCols = numel(fields);
elseif ismatrix(input)
    inputArray = input;
    [numRows, numCols] = size(inputArray);
end

output = zeros(numRows, numCols);

for i = 1:size(inputArray,2)
    minVal = min(inputArray(:,i));
    maxVal = max(inputArray(:,i));
    meanVal = mean(inputArray(:,i));
    stdVal = std(inputArray(:,i));
    Q3 = quantile(inputArray(:,i), 0.75);
    Q1 = quantile(inputArray(:,i), 0.25);
    Ql = quantile(inputArray(:,i), 0.1);
    Qr = quantile(inputArray(:,i), 0.9);
    if (Ql == 0 && Qr == 0) || (Ql == 1 && Qr == 1)
        temp = inputArray(:,i)';
        temp = mapminmax(temp, -1, +1);
        output(:,i) = temp';
    else
        switch normType
            case 'linear'
                output(:,i) = (inputArray(:,i) - minVal)/(maxVal - minVal);
            case 'mean'
                output(:,i) = (inputArray(:,i) - meanVal)/(maxVal - minVal);
            case 'std'
                output(:,i) = (inputArray(:,i) - meanVal)/stdVal;
            case 'quartile'                
                output(:,i) = (inputArray(:,i) - meanVal)/(Q3 - Q1);
            case 'quantile'                
                output(:,i) = (inputArray(:,i) - meanVal)/(Qr - Ql);
            otherwise
                error('Choose the type of normalization properly')
        end
    end
end

if isNormFirstCol == 0
    output(:,1) = inputArray(:,1);
end

if isstruct(input)
    cellArray = num2cell(output);                 
    output = cell2struct(cellArray,fields,2);
end

end