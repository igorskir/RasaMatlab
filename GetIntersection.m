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
    possibleCombs = numel(v1)*numel(v2);
    
    switch nargin
        case 2
            epsilon = defEpsilon;
        case 3
            epsilon = varargin{1,1};
        otherwise
            disp('Unexpected inputs');
            
    end
    
    % sort array elemnts for algorithms optimization
    v1 = sort(v1);
    v2 = sort(v2);

    for i = 1:numel(v1)
       for j = 1:numel(v2)           
            maxVal = max(v1(i), v2(j));         
            if ( abs(v1(i) - v2(j))/maxVal ) < epsilon
                overlap = overlap + 1;   
            else
               break; 
            end
       end
    end
    overlap = (overlap/possibleCombs)*100;
end
