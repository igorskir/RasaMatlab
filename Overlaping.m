% x = [1 , 2 , 1 , 1 , 2, 4, 4, 4, 4, 7];
% 
% y = [3, 3, 3, 2, 2 , 6, 8, 8, 8, 9, 9];
% 
% z = [1 , 5 , 1 , 1 , 5, 4, 4, 4, 4, 7];
% 
% 

function g = Overlaping(x,y,sigma)

    g = 0;
    maxM = 0;
    for i = 1 : numel(x)
        for j = 1 : numel(y)
            maxM = max(x(i),y(j));
            if ( 1 - ((abs(x(i) - y(j))) / maxM) > sigma)
                g = g + 1;
            end
        end
    end
    g = g * 100 / (numel(x)*numel(y));
end