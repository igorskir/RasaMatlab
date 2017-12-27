function outputArray = GetDataUsingModel(inputs, varargin)
    try       
        defaultRange = 'D122:U122';    
        switch nargin
            case 1
                xlRange = defaultRange;
            case 2
                xlRange = varargin{1,1};
            otherwise
                disp('Unexpected inputs');

        end
        
        % Get the model
        if ischar(varargin{1,1})
            cathDataFile = 'D:\RASA Lab\MLForCatheterDetection\!Calculations\Feature engineering.xlsm';
            sheet = 'Catheter analysis';
            model = xlsread(cathDataFile, sheet, xlRange);
        else
            model = xlRange;
        end
        numCols = numel(model);
        numRows = size(inputs, 1);
%         numRows = size(inputs, 2);
        
        if size(inputs,2) ~= numCols
            error('Input array and your model should have the same number of columns')
        end
        
        outputArray = zeros(numRows, sum(model));
        j = 1;
        for i = 1:numCols
            if model(1, i) == 1 
                outputArray(:, j) = inputs(:, i);
                j = j + 1;
            end
        end
        
    catch ME
		errorMessage = sprintf('Error in GetDataUsingModel().\n The error reported by MATLAB is:\n\n%s', ME.message);
		uiwait(warndlg(errorMessage));
		set(handles.txtInfo, 'String', errorMessage);
    end
	return;