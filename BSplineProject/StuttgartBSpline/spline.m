function spline(n,order)

% function spline(n,order)
%
% Plots the B-slpine-curve of n control-points.
% The control points can be chosen by clicking
% with the mouse on the figure.
%
% COMMAND:  spline(n,order)
% INPUT:    n     Number of Control-Points
%           order Order ob B-Splines
%                 Argnument is arbitrary
%                 default: order = 4
%
% Date:     2007-11-28
% Author:   Stefan Hüeber

close all;
if (nargin ~= 2)
	order = 4;
end
nplot = 100;

if (n < order)
	display([' !!! Error: Choose n >= order=',num2str(order),' !!!']);
	return;
end

figure(1);
hold on; box on;
set(gca,'Fontsize',16);

t = linspace(0,1,nplot);

for i = 1:n	
	title(['Choose ',num2str(i),' th. control point']);
	p(i,:) = ginput(1);
	hold off;
	plot(p(:,1),p(:,2),'k-','LineWidth',2);
	axis([0 1 0 1]);
	hold on; box on;
	if (i  >= order) 
		T = linspace(0,1,i-order+2);
		y = linspace(0,1,1000);
		p_spl = DEBOOR(T,p,y,order);
		plot(p_spl(:,1),p_spl(:,2),'b-','LineWidth',4);
	end
	plot(p(:,1),p(:,2),'ro','MarkerSize',10,'MarkerFaceColor','r');
end

title(['B-Spline-curve with ',num2str(n),' control points of order ',num2str(order)]);