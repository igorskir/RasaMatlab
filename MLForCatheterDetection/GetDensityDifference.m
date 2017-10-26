function score = GetDensityDifference(cathDataArr, tissueDataArr, varargin)
%HOW TO USE
%3 ways to call the fucntion 
% In Auto mode: 
%                  GetDensityDifference(array1, array2)
%                                   OR
%              GetDensityDifference(array1, array2, 'auto')
% In manual mode:
%              GetDensityDifference(array1, array2, 'manual', numBins)
% Example:
% areaScore = GetDensityDifference(dataCatheter(:,2), dataTissue(:,2), 'auto')
    scrSz = get(0, 'Screensize');
    try
        % Check number of inputs: manual or auto
        switch nargin
            case {2,3}
                if isempty(varargin) || strcmp(varargin{1,1}, 'auto')
                    dataArray = cat(1, cathDataArr, tissueDataArr);
                    [N, ~] = histcounts(dataArray, 'BinMethod', 'fd');
                    nBins = numel(N);
                end
            case 4
                if strcmp(varargin{1,1}, 'manual')
                    nBins = varargin{1,2};
                end
        end

        % Transpose the arrays
        if size(cathDataArr,2) > 1
            cathDataArr = cathDataArr';
        end

        if size(tissueDataArr,2) > 1
            tissueDataArr = tissueDataArr'; 
        end

        [cathBinElements, ~] = histcounts(cathDataArr, nBins);
        [tissueBinElements, ~] = histcounts(tissueDataArr, nBins);
        cathBinElementsNorm = cathBinElements/size(cathDataArr,1);
        tissueBinElementsNorm = tissueBinElements/size(tissueDataArr,1);

        score = 0;
        for i = 1:nBins
            currentDelta = cathBinElementsNorm(i) - tissueBinElementsNorm(i);
            if cathBinElements(i) ~= 0 && currentDelta > 0
                score = score + currentDelta; 
            end    
        end

        % Draw plot and tune its settings
        if nargin == 4
            mode = 'manual';
        else
            mode = 'auto';
        end
        hFig = figure;
        ax = axes('Parent', hFig);
        plot(cathBinElementsNorm, 'LineWidth', 2, 'Color', 'r');
        hold on;
        plot(tissueBinElementsNorm, 'LineWidth', 2, 'Color', 'b');
        str = sprintf('Final score in %s mode: %.4f', mode, score);
        addTitle(str)
        hold off;
        legend('Catheter PDF','Tissue PDF')
        xlabel('Data', 'FontName', 'Times New Roman')
        ylabel('PDF', 'FontName', 'Times New Roman')
        set(ax,'FontName','Times New Roman','FontSize',12);
        grid on
        set(gcf, 'Position', [1, scrSz(2), scrSz(3), scrSz(4)],...
         'Color', 'w', 'name', 'Score', 'numbertitle', 'off');

    catch ME
        errorMessage = sprintf('Error in GetDensityDifference().\n The error reported by MATLAB is:\n\n%s', ME.message);
        uiwait(warndlg(errorMessage));
        set(handles.txtInfo, 'String', errorMessage);
    end
    return;