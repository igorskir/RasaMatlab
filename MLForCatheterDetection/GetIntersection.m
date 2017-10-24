function overlap = GetIntersection(v1, v2, varargin)
    
    try
        defaultEpsilon = 0.01;
        overlap = 0;
        possibleCombs = numel(v1)*numel(v2);

        switch nargin
            case 2
                epsilon = defaultEpsilon;
            case 3
                epsilon = varargin{1,1};
            otherwise
                disp('Unexpected inputs');

        end

        % Sort array elemnts for algorithms optimization
        v1 = sort(v1);
        v2 = sort(v2);

        for i = 1:numel(v1)
           for j = 1:numel(v2)           
                maxVal = max(v1(i), v2(j));         
                if (abs(v1(i) - v2(j))/maxVal) < epsilon
                    overlap = overlap + 1;   
                elseif (v1(i) < v2(j))
                   break; 
                end
           end
        end
        overlap = (overlap/possibleCombs)*100;
    catch ME
		errorMessage = sprintf('Error in GetWeights().\n The error reported by MATLAB is:\n\n%s', ME.message);
		uiwait(warndlg(errorMessage));
		set(handles.txtInfo, 'String', errorMessage);
    end
	return;
