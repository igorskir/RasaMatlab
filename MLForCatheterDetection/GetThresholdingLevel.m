function [level, EM] = GetThresholdingLevel(img, type)

k = 1.10;

switch type
    case 'Otsu'
        [level, EM] = graythresh(img);
    case 'ContrastImproving'
        processedImg = ContrastImproving(img);
        level = graythresh(processedImg);
    case 'k*ContrastImproving'
        processedImg = ContrastImproving(img);
        level = 1.3*graythresh(processedImg);
    case 'ImprovedOtsu'
        [level, EM] = graythresh(img);
        level = k*level;
    case 'Adaptive'
        level = adaptthresh(img);
    case 'OtsuHistogram'
        [counts,~] = imhist(img,255);
        level = otsuthresh(counts);
end

end

