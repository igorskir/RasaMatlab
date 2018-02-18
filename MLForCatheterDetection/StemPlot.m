function StemPlot(edges, bins)

figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create stem
stem(edges,bins,'Marker','none','LineWidth',2,'Color',[0 0 1]);

% Create ylabel
ylabel('Elements in one bin');
ylim(axes1,[0 1000]);

% Create xlabel
xlabel('Standard deviation');
xlim(axes1,[0 60]);

box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontName','Times New Roman','FontSize',34);
