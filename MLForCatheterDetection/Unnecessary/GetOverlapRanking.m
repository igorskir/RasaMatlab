function output = GetOverlapRanking(arr1, arr2, varargin)
    try
        defEpsilon = 0.01;
               
        if size(arr1,2) ~= size(arr2,2)
            error('Inputs and Targets should have the same number of entities')
        end

        switch nargin
            case 2
                epsilon = defEpsilon;
            case 3
                epsilon = varargin{1,1};
            otherwise
                disp('Unexpected inputs');
        end
        
        numEntities = size(arr1,2);
        output = zeros(1, numEntities);
        
        for i = 1:numEntities
            output(1, i) = GetIntersection(arr1(:, i), arr2(:, i), epsilon);
        end
        
    catch ME
		errorMessage = sprintf('Error in GetWeights().\n The error reported by MATLAB is:\n\n%s', ME.message);
		uiwait(warndlg(errorMessage));
		set(handles.txtInfo, 'String', errorMessage);
    end
	return;