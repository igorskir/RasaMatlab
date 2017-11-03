function [ result] = IsDiscrete( data )
% returns true if all values of data are discrete, false otherwise
result = all(abs(floor(data)-data)==zeros(size(data)));
end

