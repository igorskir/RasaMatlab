figure1 = figure;
axes1 = axes('Parent',figure1);
boxplot(acc,'Notch','off','Labels',{'?/4', '?/8', '?/16', '?/32'}, 'MedianStyle', 'line', 'FullFactors', 'on')
xlim(axes1,[0.5 4.5]);
ylim(axes1,[0.65 1.0]);
xlabel('\Delta\phi')
ylabel('Accuracy');
box(axes1,'on');
set(axes1,'FontName','Times New Roman','FontSize', 36,'TickLabelInterpreter',...
    'none','XTick',[1 2 3 4],'XTickLabel',{'?/4','?/8','?/16','?/32'},'YGrid',...
    'on');
hold on
% plot(mean(acc(:,1)), 'dg')
plot(mean(acc(:,2)), '*',...
    'LineWidth', 1, ... 
    'MarkerEdgeColor', 'g', ...
    'MarkerSize', 10)
hold off
%%
figure1 = figure;
axes1 = axes('Parent',figure1);
boxplot(pTime2,'Notch','off','Labels',{'?/4', '?/8', '?/16', '?/32'})
xlim(axes1,[0.5 4.5]);
ylim(axes1,[0.0 15]);
xlabel('\Delta\phi')
% ylabel('Processing time');
box(axes1,'on');
set(axes1,'FontName','Times New Roman','FontSize', 36,'TickLabelInterpreter',...
    'none','XTick',[1 2 3 4],'XTickLabel',{'?/4','?/8','?/16','?/32'},'YGrid',...
    'on');