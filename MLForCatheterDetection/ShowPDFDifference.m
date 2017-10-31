function [ ] = ShowPDFDifference(dataTissue, dataCath, binCountTissue, binCountCath, isDiscrete)
scrSz = get(0, 'Screensize');
sum = 0;
[xTissue, probTissue] = GetPDFSpline(dataTissue, binCountTissue, isDiscrete);
[xCath, probCath] = GetPDFSpline(dataCath, binCountCath, isDiscrete);

tissueMin = min(xTissue);
tissueMax = max(xTissue);
cathMin = min(xCath);
cathMax = max(xCath);
z=0;
if(~isDiscrete) % non-discrete
    step = (xCath(end) - xCath(1))/(binCountCath*20);
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
else
   step=1;
   for i = xCath(1):step:xCath(end)
            cathValue=probCath(i-cathMin+1);
            if (i >= tissueMin && i <= tissueMax)
                   tissueValue=probTissue(i-tissueMin+1);
                if (tissueValue < cathValue)
                    sum = sum + cathValue-tissueValue;
                    z=i;
                end
            else
                sum = sum + cathValue;
            end
    end
end
z
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
sum

figure;
x=xCath(1):step:xCath(end);
y=spline(xCath, probCath,x);
plot(x,y);
hold on;


x=xTissue(1):step:xTissue(end);
y=spline(xTissue, probTissue,x);
plot(x,y);
legend('Catheter PDF', 'Tissue PDF')
end

