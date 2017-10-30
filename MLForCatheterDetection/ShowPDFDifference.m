function [ ] = ShowPDFDifference(dataTissue, dataCath, binCountTissue, binCountCath)
scrSz = get(0, 'Screensize');
[xTissue, probTissue] = GetPDFSpline(dataTissue,binCountTissue);
[xCath, probCath] = GetPDFSpline(dataCath,binCountCath);
step = (xCath(end) - xCath(1))/(binCountCath*20);
sum = 0;
tissueMin = min(xTissue);
tissueMax = max(xTissue);
for i = xCath(1):step:xCath(end)
    splineCath = spline(xCath, probCath,i);
    if (i >= tissueMin && i <= tissueMax)
        splineTissue = spline(xTissue, probTissue, i);
        if (splineTissue <= splineCath)
            sum = sum + (splineCath-splineTissue) * step;
        end
    else
        sum = sum + splineCath * step;
    end
end
 % Draw plot and tune its settings
hFig = figure;
ax = axes('Parent', hFig);
plot(xCath, probCath, 'LineWidth', 2, 'Color', 'b');
hold on;
plot(xTissue, probTissue, 'LineWidth', 2, 'Color', 'r');
str = sprintf('Final score: %.4f', sum);
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

