%% BiErfFit Example #1: VO2Er_transition
% 
% Copyright: Herianto Lim (http://heriantolim.com)
% Licensing: GNU General Public License v3.0
% First created: 01/11/2018
% Last modified: 02/11/2018

% Add the required packages using MatVerCon.
% addpackage('MatCommon','MatGraphics','BiErfFit');

% Clear workspace variables.
clear;

% Load data. If fails, adjust the file path supplied to the argument.
S=load('VO2Er_transition.mat');

% Perform hysteresis-curve fitting.
S.Fit=BiErfFit(S.Data,...
	'BaselineStart',.175,...
	'CenterStart',{79,[65,55]},...
	'HeightStart',{.28,[.11,.17]},...
	'WidthStart',{9.8,[5.2,28]});

%% Plotting
% Settings.
Groot.usedefault();
Groot.usedefault('latex',8,.6);
RESOLUTION=300;
AXES_SIZE=[6,6];
DATA_LINE_COLOR={'r','b'};
MODEL_LINE_COLOR={'m','g'};

% Plot data.
xLim=S.Fit.Window;
xModel=linspace(xLim(1),xLim(2),ceil(RESOLUTION/2.54*AXES_SIZE(1)));
yModel=S.Fit.model(xModel);

% Figure.
fig=docfigure(AXES_SIZE);

% Axes.
pos=[0,0,AXES_SIZE];
ax=axes('Position',pos,'XLim',xLim);
xlabel('Temperature ($^\circ$C)');
ylabel('Fresnel reflectance');
fixticklength(.2);

% Plots.
for j=1:2
	errorbar(S.Data{j}(1,:),S.Data{j}(2,:),S.Uncertainty{j},...
		'Color',DATA_LINE_COLOR{j},'LineWidth',.5);
end
for j=1:2
	plot(xModel,yModel{j},'Color',MODEL_LINE_COLOR{j},'LineWidth',.7);
end

% Reconfigure the layout.
margin=ax.TightInset+.1;
ax.Position=pos+[margin(1),margin(2),0,0];
pos=pos+[0,0,margin(1)+margin(3),margin(2)+margin(4)];
set(fig,{'Position','PaperPosition','PaperSize'},{pos,pos,pos(3:4)});

% Saving.
print(fig,'VO2Er_transition.png','-dpng',sprintf('-r%d',RESOLUTION));
close(fig);
