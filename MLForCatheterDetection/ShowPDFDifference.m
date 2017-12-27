function [ area ] = ShowPDFDifference(dataTissue, dataCath, binCountTissue, binCountCath, isDiscrete, varargin)

if nargin == 5 || nargin == 6
    scrSz = get(0, 'Screensize');
    area = 0;
    [xTissue, probTissue] = GetPDFSpline(dataTissue, binCountTissue, isDiscrete);
    [xCath, probCath] = GetPDFSpline(dataCath, binCountCath, isDiscrete);

    tissueMin = min(xTissue);
    tissueMax = max(xTissue);
    cathMin = min(xCath);
    cathMax = max(xCath);
    splineMethod = 'pchip';
    if(~isDiscrete) % non-discrete
        step = (xCath(end) - xCath(1))/(binCountCath*20);
        for i = xCath(1):step:xCath(end)
            %splineCath = spline(xCath, probCath,i);
            splineCath = interp1(xCath, probCath,i,splineMethod);
            if (i >= tissueMin && i <= tissueMax)
                %splineTissue = spline(xTissue, probTissue, i);
                splineTissue = interp1(xTissue, probTissue, i,splineMethod);
                if (splineTissue <= splineCath)
                    area = area + (splineCath-splineTissue) * step;
                end
            else
                area = area + splineCath * step;
            end
        end
        maxScore=sum(interp1(xCath, probCath,xCath(1):step:xCath(end),splineMethod))*step;
        area=area/maxScore;
    else
       step=1;
       for i = xCath(1):step:xCath(end)
                cathValue=probCath(i-cathMin+1);
                if (i >= tissueMin && i <= tissueMax)
                       tissueValue=probTissue(i-tissueMin+1);
                    if (tissueValue < cathValue)
                        area = area + cathValue-tissueValue;
                    end
                else
                    area = area + cathValue;
                end
       end
    end
end

if nargin == 6
    % Draw plot and tune its settings
    if ~isempty(varargin{1,1}) && strcmp(varargin{1,1}, 'uPlot')
        hFig = figure;
        ax = axes('Parent', hFig);
        plot(xCath, probCath, 'LineWidth', 2, 'Color', 'b');
        hold on;
        plot(xTissue, probTissue, 'LineWidth', 2, 'Color', 'r');
    end

    % Smoothed plot
    if ~isempty(varargin{1,1}) && strcmp(varargin{1,1}, 'sPlot')
        hFig = figure;
        ax = axes('Parent', hFig);
        x = xCath(1):step:xCath(end);
        y = spline(xCath, probCath, x);
        y = interp1(xCath, probCath, x, splineMethod);
        plot(x, y, 'LineWidth', 2, 'Color', 'b');
        hold on;
        x=xTissue(1):step:xTissue(end);
        y=interp1(xTissue, probTissue, x, splineMethod);
        plot(x,y, 'LineWidth', 2, 'Color', 'r');
    end
    
    str = sprintf('Final score: %.4f', area);
    AddTitle(str)
    hold off;
    legend('Catheter PDF', 'Tissue PDF')
    xlabel('Data', 'FontName', 'Times New Roman');
    ylabel('PDF', 'FontName', 'Times New Roman');
    set(ax,'FontName','Times New Roman','FontSize',12);
    grid on
    set(gcf, 'Position', [1, scrSz(2), scrSz(3), scrSz(4)],...
    'Color', 'w', 'name', 'Score', 'numbertitle', 'off');
end
end
