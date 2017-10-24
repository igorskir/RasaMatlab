% x = [1, 2, 1, 1, 2, 4, 4, 4, 4, 7];
% 
% y = [3, 3, 3, 2, 2, 6, 8, 8, 8, 9, 9];
% 
% z = [1, 5, 1, 1, 5, 4, 4, 4, 4, 7];
% 
% 

function overlap = GetIntersection(v1, v2, varargin)
    
    defEpsilon = 0.01;
    overlap = 0;
    
    switch nargin
        case 2
            epsilon = defEpsilon;
        case 3
            epsilon = varargin{1,1};
        otherwise
            disp('Unexpected inputs');
            
    end

    for i = 1:numel(v1)
       for j = 1:numel(v2)           
            minVal = min(v1(i), v2(j));
            maxVal = max(v1(i), v2(j));
            if (1 - (minVal/maxVal)) < epsilon
                overlap = overlap + 1;            
            end
       end
    end
    overlap = (overlap/(numel(v1)*numel(v2)))*100;
end
