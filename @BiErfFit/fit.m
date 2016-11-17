function obj=fit(obj)
%% BiErf Fitting
%  This method performs a curve fitting to a hysteresis data with a linear
%  combination of error functions using the MATLAB Curve Fitting toolbox.
%
%  The fit results and statistics will be stored in the object properties after
%  a successful calling of this method.
%
%  The object properties related to the start points and constraints for the fit
%  parameters will also be initialized.
%
%  This algorithm is far from perfect. A lot of things need fixing to adapt to
%  different scenarios or be more efficient.
%
% Requires package:
%  - PeakFit_v1.0.0+
%
% Tested on:
%  - MATLAB R2015b
%
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 25/05/2013
% Last modified: 17/11/2016

x=obj.XData;
y=obj.YData;
if isempty(x) || isempty(y)
	return
end

%% Defaults
K=-1000;
BASELINE_TOL=.225;

%% Initialize the Fit Data
for j=1:2
	assert(numel(x{j})==numel(y{j}),...
		'BiErfFit:init:InconsistentNumPoints',...
		'The number of X and Y data points must be equal.');
end

% Sort the points (x,y) in ascending order
for j=1:2
	[x{j},ix]=sort(x{j});
	y{j}=y{j}(ix);
end

% Update the object data
obj.XData=x;
obj.YData=y;

% Trim data points
if ~isempty(obj.Window)
	for j=1:2
		ix=x{j}>=obj.Window(1) & x{j}<=obj.Window(2);
		assert(sum(ix)>=obj.MinNumPoints,...
			'BiErfFit:init:InsufficientNumPoints',...
			['The number of data points within the fit window must be at ',...
				'least %d.'],obj.MinNumPoints);
		x{j}=x{j}(ix);
		y{j}=y{j}(ix);
	end
end

% Perform linear mapping from [x(1),x(end)]->[0,1] and [min(y),max(y)]->[0,1]
n=[numel(x{1}),numel(x{2})];
c1=min(x{1}(1),x{2}(1));
c2=max(x{1}(n(1)),x{2}(n(2)));
h1=min(min(y{1}),min(y{2}));
h2=max(max(y{1}),max(y{2}));
a=1/(c2-c1);
b=-a*c1;
c=1/(h2-h1);
d=-c*h1;
for j=1:2
	x{j}=a*x{j}+b;
	y{j}=c*y{j}+d;
end

% Update the fit window
if isempty(obj.Window)
	obj.Window=[c1,c2];
end

%% Initialize the Fit Parameters
ne=obj.NumErfs;
center=cell(3,2);
height=cell(3,2);
width=cell(3,2);

if ~isempty(obj.CenterStart)
	center(1,:)=obj.CenterStart;
end
if ~isempty(obj.CenterLow)
	center(2,:)=obj.CenterLow;
end
if ~isempty(obj.CenterUp)
	center(3,:)=obj.CenterUp;
end
if ~isempty(obj.HeightStart)
	height(1,:)=obj.HeightStart;
end
if ~isempty(obj.HeightLow)
	height(2,:)=obj.HeightLow;
end
if ~isempty(obj.HeightUp)
	height(3,:)=obj.HeightUp;
end
if ~isempty(obj.WidthStart)
	width(1,:)=obj.WidthStart;
end
if ~isempty(obj.WidthLow)
	width(2,:)=obj.WidthLow;
end
if ~isempty(obj.WidthUp)
	width(3,:)=obj.WidthUp;
end

% Find the number of erf functions
for i=1:3
	for j=1:2
		ne(j)=max([ne(j),...
			numel(center{i,j}),numel(height{i,j}),numel(width{i,j})]);
	end
end

% Transform the constraints, and fill the blanks with NaN
for i=1:3
	for j=1:2
		center{i,j}=[a*center{i,j}+b,nan(1,ne(j)-numel(center{i,j}))];
		height{i,j}=[c*height{i,j},nan(1,ne(j)-numel(height{i,j}))];
		width{i,j}=[a*width{i,j},nan(1,ne(j)-numel(width{i,j}))];
	end
end

% If lower bounds exceed upper bounds, swap them
for j=1:2
	for i=1:ne(j)
		if center{2,j}(i)>center{3,j}(i)
			c2=center{2,j}(i);
			center{2,j}(i)=center{3,j}(i);
			center{3,j}(i)=c2;
		end
		if height{2,j}(i)>height{3,j}(i)
			h2=height{2,j}(i);
			height{2,j}(i)=height{3,j}(i);
			height{3,j}(i)=h2;
		end
		if width{2,j}(i)>width{3,j}(i)
			w2=width{2,j}(i);
			width{2,j}(i)=width{3,j}(i);
			width{3,j}(i)=w2;
		end
	end
end

% Fill the blank constraints with estimates from peak fitting of the gradient.
% This block needs revision, as the derivative may contain NaN values or spikes.
% Suggestions: spline interpolation, derivative approximations with higher order
% terms, spike removal with statistical filter.
for j=1:2
	x1=mean([x{j}(1:n(j)-1);x{j}(2:n(j))]);
	y1=diff(y{j})./diff(x{j});
	h1=1-2*(abs(max(y1)/min(y1))<1);% +1 if peak is upwards, -1 if downwards
	
	try
		Fit=PeakFit(x1,h1*y1,...
			'NumPeaks',ne(j),...
			'PeakShape','Gaussian',...
			'CenterStart',center{1,j},...
			'CenterLow',center{2,j},...
			'CenterUp',center{3,j},...
			'AreaStart',h1*height{1,j},...
			'AreaLow',h1*height{2,j},...
			'AreaUp',h1*height{3,j},...
			'WidthStart',width{1,j},...
			'WidthLow',width{2,j},...
			'WidthUp',width{3,j},...
			'BaselinePolyOrder',-1);
	catch ME1
		if isempty(regexpi(ME1.identifier,':InvalidInput$','once'))
			continue
		else
			rethrow(ME1);
		end
	end
	
	for i=1:3
		ix=~isfinite(center{i,j}) & isfinite(Fit.Center(i,:));
		center{i,j}(ix)=Fit.Center(i,ix);
		ix=~isfinite(height{i,j}) & isfinite(Fit.Area(i,:));
		height{i,j}(ix)=h1*Fit.Area(i,ix);
		ix=~isfinite(width{i,j}) & isfinite(Fit.Width(i,:));
		width{i,j}(ix)=Fit.Width(i,ix);
	end
end

% Fill the blank start points and constraints for the center, height, and width
w3=obj.MovMeanWidth;
if w3>=1
	w3=w3/(n(1)+n(2));
else
	w3=w3/2;
end
h1=zeros(2);
h2=2*obj.TolFun;
for j=1:2
	h1(1,j)=mean(y{j}(x{j}<=w3));
	h1(2,j)=mean(y{j}(x{j}>=1-w3));
	ix=(y{j}>h1(1,j)+BASELINE_TOL & y{j}<h1(2,j)-BASELINE_TOL) ...
		| (y{j}>h1(2,j)+BASELINE_TOL & y{j}<h1(1,j)-BASELINE_TOL);
	w2=mean(diff(x{j}));
	w1=sum(ix)*w2;
	c1=x{j}(round(mean(find(ix))));
	
	center{1,j}(~isfinite(center{1,j}))=c1;
	center{2,j}(~isfinite(center{2,j}))=c1-w1;
	center{3,j}(~isfinite(center{3,j}))=c1+w1;
	
	height{1,j}(~isfinite(height{1,j}))=(h1(2,j)-h1(1,j))/sqrt(ne(j));
	for i=1:ne(j)
		if height{1,j}(i)<0
			if ~isfinite(height{2,j}(i))
				height{2,j}(i)=-1;
			end
			if ~isfinite(height{3,j}(i))
				height{3,j}(i)=-h2;
			end
		else
			if ~isfinite(height{2,j}(i))
				height{2,j}(i)=h2;
			end
			if ~isfinite(height{3,j}(i))
				height{3,j}(i)=1;
			end
		end
	end
	
	width{1,j}(~isfinite(width{1,j}))=w1;
	width{2,j}(~isfinite(width{2,j}))=w2;
	width{3,j}(~isfinite(width{3,j}))=2*w1;
end

% Fill the blank start points and constraints for the baseline
h1=mean(h1(1,:));
baseline=h1+[0;-BASELINE_TOL;BASELINE_TOL];
if ~isempty(obj.BaselineStart)
	baseline(1)=c*obj.BaselineStart+d;
end
if ~isempty(obj.BaselineLow)
	baseline(2)=c*obj.BaselineLow+d;
end
if ~isempty(obj.BaselineUp)
	baseline(3)=c*obj.BaselineUp+d;
end

% Ensure that the start points and constraints for the height satisfy the
% boundary conditions: sum(height{:,1})=sum(height{:,2})
h3=[mean([sum(height{1,1}),sum(height{1,2})]);
	min(sum(height{2,1}),sum(height{2,2}));
	max(sum(height{3,1}),sum(height{3,2}))];
if h3(1)<0
	h3(3)=min(h3(3),-h2);
else
	h3(2)=max(h3(2),h2);
end
for i=1:3
	for j=1:2
		height{i,j}(1)=h3(i)-sum(height{i,j}(2:ne(j)));
	end
end

%% Construct the Fit Coefficients and Expression
problem='k';

coeff=[{'C','H','cA1','cB1','wA1','wB1'},cell(1,3*(sum(ne)-2))];
for i=2:ne(1)
	coeff{5+i}=sprintf('hA%d',i);
	coeff{3+ne(1)+ne(2)+i}=sprintf('cA%d',i);
	coeff{1+2*ne(1)+2*ne(2)+i}=sprintf('wA%d',i);
end
for i=2:ne(2)
	coeff{4+ne(1)+i}=sprintf('hB%d',i);
	coeff{2+2*ne(1)+ne(2)+i}=sprintf('cB%d',i);
	coeff{3*ne(1)+2*ne(2)+i}=sprintf('wB%d',i);
end

expr=sprintf('%s+%s',coeff{1:2});
if ne(1)==1
	expr=[expr,sprintf('+%s/2*(1+erf(2*sqrt(log(2))*(x-%s)/%s))',...
		coeff{2},coeff{3},coeff{5})];
else
	expr=[expr,'+(',coeff{2},sprintf('-%s',coeff{7:5+ne(1)}),')',...
		sprintf('/2*(1+erf(2*sqrt(log(2))*(x-%s)/%s))',...
			coeff{3},coeff{5})];
	for i=2:ne(1)
		expr=[expr,sprintf('+%s/2*(1+erf(2*sqrt(log(2))*(x-%s)/%s))',...
			coeff{5+i},coeff{3+ne(1)+ne(2)+i},...
			coeff{1+2*ne(1)+2*ne(2)+i})]; %#ok<AGROW>
	end
end
if ne(2)==1
	expr=[expr,sprintf('-%s/2*(1+erf(2*sqrt(log(2))*(x-%s+%s)/%s))',...
		coeff{2},problem,coeff{4},coeff{6})];
else
	expr=[expr,'-(',coeff{2},sprintf('-%s',coeff{6+ne(1):4+ne(1)+ne(2)}),')',...
		sprintf('/2*(1+erf(2*sqrt(log(2))*(x-%s+%s)/%s))',...
			problem,coeff{4},coeff{6})];
	for i=2:ne(2)
		expr=[expr,sprintf('-%s/2*(1+erf(2*sqrt(log(2))*(x-%s+%s)/%s))',...
			coeff{4+ne(1)+i},problem,coeff{2+2*ne(1)+ne(2)+i},...
			coeff{3*ne(1)+2*ne(2)+i})]; %#ok<AGROW>
	end
end

%% Further Transformation before Fitting
% Reflect the second curve about an axis x=K/2
x{2}=K-x{2};
[x{2},ix]=sort(x{2});
y{2}=y{2}(ix);

% Combine the two curves together
x=[x{1},x{2}];
y=[y{1},y{2}];

%% MATLAB Curve Fitting ToolBox
FitConstraint=cell(3,1);
for i=1:3
	FitConstraint{i}=[baseline(i),sum(height{i,1}),...
		center{i,1}(1),center{i,2}(1),width{i,1}(1),width{i,2}(1),...
		height{i,1}(2:ne(1)),height{i,2}(2:ne(2)),...
		center{i,1}(2:ne(1)),center{i,2}(2:ne(2)),...
		width{i,1}(2:ne(1)),width{i,2}(2:ne(2))];
end
FitOption=fitoptions('Method',obj.Method,...
	'Robust',obj.Robust,...
	'StartPoint',FitConstraint{1},...
	'Lower',FitConstraint{2},...
	'Upper',FitConstraint{3},...
	'Algorithm',obj.Algorithm,...
	'DiffMaxChange',obj.DiffMaxChange,...
	'DiffMinChange',obj.DiffMinChange,...
	'MaxFunEvals',obj.MaxFunEvals,...
	'MaxIter',obj.MaxIters,...
	'TolFun',obj.TolFun,...
	'TolX',obj.TolX,...
	'Display','off');
FitType=fittype(expr,...
	'coefficients',coeff,...
	'problem',problem,...
	'options',FitOption);
warning('off','all');
[FitObj,gof,FitOutput]=fit(x',y',FitType,'problem',K);
warning('on','all');
fitResult=[coeffvalues(FitObj);confint(FitObj)];

%% Update Object Properties
% Update the number of error functions
obj.NumErfs=ne;

% Reverse mapping coefficients
a=1/a;
b=-a*b;
c=1/c;
d=-c*d;

% Update the constraint parameters, after a reverse mapping
for i=1:3
	for j=1:2
		center{i,j}=a*center{i,j}+b;
		height{i,j}=c*height{i,j};
		width{i,j}=a*width{i,j};
	end
end
baseline=c*baseline+d;
obj.CenterStart=center(1,:);
obj.CenterLow=center(2,:);
obj.CenterUp=center(3,:);
obj.HeightStart=height(1,:);
obj.HeightLow=height(2,:);
obj.HeightUp=height(3,:);
obj.WidthStart=width(1,:);
obj.WidthLow=width(2,:);
obj.WidthUp=width(3,:);
obj.BaselineStart=baseline(1);
obj.BaselineLow=baseline(2);
obj.BaselineUp=baseline(3);

% Update the fit output, after a reverse mapping
h1=c*fitResult(:,2);
c1=a*fitResult(:,3:4)+b;
w1=a*fitResult(:,5:6);
h2=c*fitResult(:,7:5+ne(1));
h3=c*fitResult(:,6+ne(1):4+ne(1)+ne(2));
c2=a*fitResult(:,5+ne(1)+ne(2):3+2*ne(1)+ne(2))+b;
c3=a*fitResult(:,4+2*ne(1)+ne(2):2+2*ne(1)+2*ne(2))+b;
w2=a*fitResult(:,3+2*ne(1)+2*ne(2):1+3*ne(1)+2*ne(2));
w3=a*fitResult(:,2+3*ne(1)+2*ne(2):3*ne(1)+3*ne(2));
obj.Baseline=c*fitResult(:,1)+d;
obj.Height={[h1-sum(h2,2),h2],[h1-sum(h3,2),h3]};
obj.Center={[c1(:,1),c2],[c1(:,2),c3]};
obj.Width={[w1(:,1),w2],[w1(:,2),w3]};

% Update the goodness of fit
obj.RelStDev=gof.rmse;% due to the transformation to [0,1],
                      % the rmse is effectively a relative standard deviation
obj.CoeffDeterm=gof.rsquare;
obj.AdjCoeffDeterm=gof.adjrsquare;

% Update the performance stats of the fitting
obj.NumFunEvals=FitOutput.funcCount;
obj.NumIters=FitOutput.iterations;
obj.ExitFlag=FitOutput.exitflag;

end