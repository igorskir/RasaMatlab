function DrawHists(v1, v2, varargin)

    switch nargin
        case 2
            numBins1 = numel(v1);
            numBins2 = numel(v2);
        case 3
            numBins1 = varargin{1,1};
            numBins2 = varargin{1,1};
        otherwise
            disp('Unexpected inputs');
    end


    histogram(v2, numBins2, 'FaceAlpha', 0.5);
    hold on;
    histogram(v1, numBins1, 'FaceAlpha', 0.5);
    grid on;
    histogram(v2, 11987, 'FaceAlpha', 0.5);
    hold on;
    histogram(v1, 881, 'FaceAlpha', 0.5);
    grid on;
    % Put up legend.
    legend1 = sprintf('Mean = %.1f', mean(v1));
    legend2 = sprintf('Mean = %.1f', mean(v2));
    legend({legend1, legend2});
end