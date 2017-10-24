function [out, fn] = struct2mat(input)

fn = fieldnames(input);
numCols = numel(fn);
numRows = numel(input);
out = zeros(numRows, numCols);

for i = 1:numCols  
    for j = 1:numRows         
%         out(j, i) = getfield(input(j), fn{i});
        out(j, i) = input(j).(fn{i});
    end
end