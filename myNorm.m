function arr = myNorm(x)
    arr = zeros(1, numel(x));
    stdLocal = std(x);
    meanLocal = mean(x);
    for i = 1 : numel(x)
        arr(i) = (x(i) - meanLocal)/stdLocal;
    end
end