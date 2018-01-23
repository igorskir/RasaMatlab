numIterations = 150;
processingTime = zeros(numIterations, 1);
str = 1;
for i = 1:numIterations
    tic;
    reconstruct3Dcath
    processingTime(str) = toc;
    str = str + 1;
end

