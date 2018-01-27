function [yModel,yErf,yBaseline]=model(obj,varargin)
%% Create Fit Model Using the Fit Results
%  [yModel,yPeak,yBaseline]=obj.model() returns a set of fit models evaluated at
%  the data points.
%
%  [yModel,yPeak,yBaseline]=obj.model(x) evaluates the fit models at x instead.
%
%  [yModel,yPeak,yBaseline]=obj.model(x_?,x_?) evaluates the fit models for
%  curve ? at x_? and curve ? at x_?.
%
% Outputs:
%  yModel: The reconstructed data points from the fit results.
%
%  yPeak: A cell containing the resconstructed data points of each peak.
%
%  yBaseline: The reconstructed baseline points.
%
% Requires package:
%  - Common_v1.0.0+
%
% Tested on:
%  - MATLAB R2013b
%  - MATLAB R2015b
%
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 29/05/2013
% Last modified: 25/10/2016

yModel=cell(1,2);
yErf=cell(1,2);
yBaseline=cell(1,2);

center=obj.Center;
height=obj.Height;
width=obj.Width;
baseline=obj.Baseline;

if isempty(center) || isempty(height) || isempty(width) || isempty(baseline)
	return
end

%% Parse Inputs
if nargin==1
	xModel=obj.XData;
elseif nargin>2
	error('BiErfFit:model:TooManyInput',...
		'At most one input argument is accepted.');
elseif isrealvector(varargin{1})
	xModel=repmat(varargin(1),1,2);
elseif iscell(varargin{1}) && numel(varargin{1})==2 ...
		&& isrealvector(varargin{1}{1}) && isrealvector(varargin{1}{2})
	xModel=varargin{1};
else
	error('BiErfFit:model:InvalidInput',...
		['Input to the domains for the model must be either a real vector ',...
			'or a cell containing two real vectors.']);
end

%% Construct Baseline
for i=1:2
	yBaseline{i}=baseline(1,1)*ones(size(xModel{i}));
end

%% Construct Erf
N=obj.NumErfs;
for i=1:2
	yErf{i}=cell(1,N(i));
	for j=1:N(i)
		yErf{i}{j}=BiErfFit.fnerf(xModel{i},...
			center{i}(1,j),height{i}(1,j),width{i}(1,j));
	end
end

%% Construct Model
for i=1:2
	yModel{i}=yBaseline{i};
	for j=1:N(i)
		yModel{i}=yModel{i}+yErf{i}{j};
	end
end

end