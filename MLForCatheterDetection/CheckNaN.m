function numNaN = CheckNaN(inp)

if isstruct(inp)
    inp = struct2mat(inp);
end

[numRows, numCols] = size(inp);
numNaN = 0;
for i = 1:numRows
    for j = 1:numCols
       if isnan(inp(i,j))
        numNaN = numNaN + 1;
       end
    end
end
end