function [probDistCatheter, probDistTissue] = CreateProbabilityPlot(catheterData,tissueData, featureName)

% Force all inputs to be column vectors
scrSz = get(0, 'Screensize');
catheterData = catheterData(:);
tissueData = tissueData(:);

% Prepare figure
clf;
hold on;
% tab1 = uitab('Title','Tab1');
% ax1 = axes(tab1);
LegHandles = []; LegText = {};
probplot('normal'); % create empty plot of desired type
title(featureName);

% Plot data originally in dataset "Area (catheter)"
hLine = probplot(gca,catheterData,[],[],'noref'); % add data to existing plot
set(hLine,'Color',[0.333333 0 0.666667],'Marker','o', 'MarkerSize',6);
xlabel('Data', 'FontName', 'Times New Roman','FontSize',14);
ylabel('Probability', 'FontName', 'Times New Roman','FontSize',14);
LegHandles(end+1) = hLine;
LegText{end+1} = strcat(featureName, ' (catheter)');

% Plot data originally in dataset "Area (tissue)"
hLine = probplot(gca,tissueData,[],[],'noref'); % add data to existing plot
set(hLine,'Color',[0.333333 0.666667 0],'Marker','o', 'MarkerSize',6);
xlabel('Data', 'FontName', 'Times New Roman','FontSize',14);
ylabel('Probability', 'FontName', 'Times New Roman','FontSize',14);
LegHandles(end+1) = hLine;
LegText{end+1} = strcat(featureName, ' (tissue)');

% Create fit "featureName (catheter)"
% This fit does not appear on the plot
probDistCatheter = fitdist(catheterData,'kernel','kernel','normal','support','unbounded');

% Create fit "featureName (tissue)"
% This fit does not appear on the plot
probDistTissue = fitdist(tissueData,'kernel','kernel','normal','support','unbounded');

% Adjust figure
box on;
grid on;
hold off;

% Create legend from accumulated handles and labels
hLegend = legend(LegHandles,LegText,'Orientation', 'vertical', 'FontName', 'Times New Roman', 'FontSize', 12, 'Location', 'southeast');
title(hLegend, 'Feature');
set(hLegend,'Interpreter','none');
set(gcf, 'Position', [scrSz(3), 0, scrSz(3), scrSz(4)], 'Color', 'w', 'name', featureName, 'numbertitle', 'off');

% -----------------------------------------------
function f = cdfnp(x,y,cens,freq,support,kernel,width)
%CDFNP Compute cdf for non-parametric fit, used in probability plot

f = ksdensity(y,x,'cens',cens,'weight',freq,'function','cdf',...
    'support',support,'kernel',kernel,'width',width);
