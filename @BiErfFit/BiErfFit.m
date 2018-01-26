classdef BiErfFit
%% BiErfFit Class
%  Fit a hysteresis curve simultaneously with two linear combinations of error
%  functions.
%
% Requires package:
%  - Common_v1.0.0+
%  - PeakFit_v1.0.0+
%
% Tested on:
%  - MATLAB R2015b
%
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 24/05/2013
% Last modified: 04/11/2016

%% Properties
% Data, start points, and constraints for fitting
properties (Dependent=true)
	Data={};
end
properties
	XData={};
	YData={};
	Window
	NumErfs=[1,1];
	CenterStart={};
	CenterLow={};
	CenterUp={};
	HeightStart={};
	HeightLow={};
	HeightUp={};
	WidthStart={};
	WidthLow={};
	WidthUp={};
	BaselineStart
	BaselineLow
	BaselineUp
end

% Fitting options
properties (Constant=true)
	Method='NonlinearLeastSquares';
end
properties
	Robust='off';
	Algorithm='Trust-Region';
	MovMeanWidth=.02;
	DiffMaxChange=.1;
	DiffMinChange=1e-8;
	MaxFunEvals=1e5;
	MaxIters=1e3;
	TolFun=1e-6;
	TolX=1e-6;
end

% Fit results
properties (SetAccess=protected)
	Center={};
	Height={};
	Width={};
	Baseline
end
properties (Dependent=true)
	Asymptote
	DynamicRange
	HysteresisWidth
	TransitionWidth
end

% Fit errors and performance
properties (SetAccess=protected)
	RelStDev=0;
	CoeffDeterm=0;
	AdjCoeffDeterm=0;
	NumFunEvals=0;
	NumIters=0;
	ExitFlag
end

% Class dictionary
properties (Constant=true,GetAccess=protected)
	MinNumPoints=10;
	RobustList={'on','off','LAR','Bisquare'};
	AlgorithmList={'Levenberg-Marquardt','Trust-Region'};
end

%% Methods
methods
	% Constructor
	function obj=BiErfFit(varargin)
		N=nargin;
		if N==0
			return
		end

		% Copy the BiErfFit object given in the first argument if any
		k=1;
		if isa(varargin{k},class(obj))
			obj=varargin{k};
			k=k+1;
		end

		% Parse input to Data points
		if k<=N && iscell(varargin{k}) && numel(varargin{k})==2 ...
				&& isrealmatrix(varargin{k}{1})...
				&& isrealmatrix(varargin{k}{2})
			if isvector(varargin{k}{1}) && isvector(varargin{k}{2})
				if k<N && iscell(varargin{k+1}) && numel(varargin{k+1})==2 ...
						&& isrealvector(varargin{k+1}{1})...
						&& isrealvector(varargin{k+1}{2})
					obj.XData=varargin{k};
					obj.YData=varargin{k+1};
					k=k+2;
				end
			else
				obj.Data=varargin{k};
				k=k+1;
			end
		end

		% Parse inputs to the parameters
		P=properties(obj);
		while k<N
			if ~isstringscalar(varargin{k})
				break
			end
			ix=strcmpi(varargin{k},P);
			if any(ix)
				obj.(P{ix})=varargin{k+1};
				k=k+2;
			else
				break
			end
		end
		assert(k>N,...
			'BiErfFit:UnexpectedInput',...
			'One or more inputs are not recognized.');

		% Perform erf fitting
		obj=fit(obj);
	end

	% Get Methods
	function d=get.Data(obj)
		x=obj.XData;
		y=obj.YData;
		if isempty(x) || isempty(y)
			d=[];
		elseif numel(x{1})~=numel(y{1}) || numel(x{2})~=numel(y{2})
			d=[];
		else
			d={[x{1};y{1}],[x{2};y{2}]};
		end
	end

	function x=get.Asymptote(obj)
		baseline=obj.Baseline;
		height=obj.Height;
		if isempty(baseline) || isempty(height)
			x=[];
		else
			x=[baseline,baseline([1,3,2])+sum(height{1},2)];
		end
	end

	function x=get.DynamicRange(obj)
		baseline=obj.Baseline;
		height=obj.Height;
		if isempty(baseline) || isempty(height)
			x=[];
		else
			x=sum(height{1},2)/baseline(1,1);
		end
	end

	function x=get.HysteresisWidth(obj)
		center=obj.Center;
		if isempty(center)
			x=[];
		elseif center{1}(1,1)<center{2}(1,1)
			x=center{2}(:,1)-center{1}([1,3,2],1);
		else
			x=center{1}(:,1)-center{2}([1,3,2],1);
		end
	end

	function x=get.TransitionWidth(obj)
		width=obj.Width;
		if isempty(width)
			x=[];
		else
			x=cell(1,2);
			for j=1:2
				[~,k]=max(width{j}(1,:));
				x{j}=width{j}(:,k);
			end
		end
	end

	% Set methods
	function obj=set.Data(obj,d)
		if isempty(d)
			obj.XData={};
			obj.YData={};
			return
		end
		ME=MException('PeakFit:InvalidInput',...
			['Input to set the Data must be a cell with two elements, both of',...
				'which must be a real matrix of size [n,2] or [2,n].']);
		if iscell(d) && numel(d)==2
			x=cell(1,2);
			y=cell(1,2);
			for i=1:2
				if isrealmatrix(d{i})
					[m,n]=size(d{i});
					if m==2
						x{i}=d{i}(1,:);
						y{i}=d{i}(2,:);
					elseif n==2
						x{i}=d{i}(:,1);
						y{i}=d{i}(:,2);
					end
				else
					throw(ME);
				end
			end
			obj.XData=x;
			obj.YData=y;
		else
			throw(ME);
		end
	end

	function obj=set.XData(obj,x)
		if isempty(x)
			obj.XData={};
		elseif iscell(x) && numel(x)==2 ...
				&& isrealvector(x{1}) && all(isfinite(x{1}))...
				&& isrealvector(x{2}) && all(isfinite(x{2}))
			if numel(x{1})<obj.MinNumPoints || numel(x{2})<obj.MinNumPoints
				error('BiErfFit:setXData:InsufficientNumPoints',...
					['Each vector in the cell input to XData must have at least ',...
						'%d elements.'],obj.MinNumPoints);
			else
				obj.XData={x{1}(:).',x{2}(:).'};
			end
		else
			error('BiErfFit:setXData:InvalidInput',...
				['Input to set the XData must be a cell containing two ',...
					'vectors of finite real numbers.']);
		end
	end

	function obj=set.YData(obj,x)
		if isempty(x)
			obj.YData={};
		elseif iscell(x) && numel(x)==2 ...
				&& isrealvector(x{1}) && all(isfinite(x{1}))...
				&& isrealvector(x{2}) && all(isfinite(x{2}))
			if numel(x{1})<obj.MinNumPoints || numel(x{2})<obj.MinNumPoints
				error('BiErfFit:setYData:InsufficientNumPoints',...
					['Each vector in the cell input to YData must have at least ',...
						'%d elements.'],obj.MinNumPoints);
			else
				obj.YData={x{1}(:).',x{2}(:).'};
			end
		else
			error('BiErfFit:setYData:InvalidInput',...
				['Input to set the YData must be a cell containing two ',...
					'vectors of finite real numbers.']);
		end
	end

	function obj=set.Window(obj,x)
		if isempty(x)
			obj.Window=[];
		elseif isrealvector(x) && numel(x)==2
			obj.Window=sort(x(:)).';
		else
			error('BiErfFit:setWindow:InvalidInput',...
				'Input to set the Window must be a real vector of length two.');
		end
	end

	function obj=set.NumErfs(obj,x)
		if isintegervector(x) && numel(x)==2 && all(x>0)
			obj.NumErfs=x(:).';
		else
			error('BiErfFit:setNumErfs:InvalidInput',...
				['Input to set the NumErfs must be a positive integer vector ',...
					'of length two.']);
		end
	end

	function obj=set.CenterStart(obj,x)
		if isempty(x)
			obj.CenterStart={};
		elseif isrealvector(x) && numel(x)==2
				obj.CenterStart={x(1),x(2)};
		elseif iscell(x) && numel(x)==2 ...
				&& isrealvector(x{1}) && isrealvector(x{2})
			obj.CenterStart={x{1}(:).',x{2}(:).'};
		else
			error('BiErfFit:setCenterStart:InvalidInput',...
				['Input to set the CenterStart must be either a two-element ',...
					'vector or a cell containing two real vectors.']);
		end
	end

	function obj=set.CenterLow(obj,x)
		if isempty(x)
			obj.CenterLow={};
		elseif isrealvector(x) && numel(x)==2
				obj.CenterLow={x(1),x(2)};
		elseif iscell(x) && numel(x)==2 ...
				&& isrealvector(x{1}) && isrealvector(x{2})
			obj.CenterLow={x{1}(:).',x{2}(:).'};
		else
			error('BiErfFit:setCenterLow:InvalidInput',...
				['Input to set the CenterLow must be either a two-element ',...
					'vector or a cell containing two real vectors.']);
		end
	end

	function obj=set.CenterUp(obj,x)
		if isempty(x)
			obj.CenterUp={};
		elseif isrealvector(x) && numel(x)==2
				obj.CenterUp={x(1),x(2)};
		elseif iscell(x) && numel(x)==2 ...
				&& isrealvector(x{1}) && isrealvector(x{2})
			obj.CenterUp={x{1}(:).',x{2}(:).'};
		else
			error('BiErfFit:setCenterUp:InvalidInput',...
				['Input to set the CenterUp must be either a two-element ',...
					'vector or a cell containing two real vectors.']);
		end
	end

	function obj=set.HeightStart(obj,x)
		if isempty(x)
			obj.HeightStart={};
		elseif isrealscalar(x)
			obj.HeightStart={x,x};
		elseif iscell(x)
			N=numel(x);
			if N==1 && isrealscalar(x{1})
				obj.HeightStart=[x(1),x(1)];
			elseif N==2 && isrealvector(x{1}) && isrealvector(x{2})
				obj.HeightStart={x{1}(:).',x{2}(:).'};
			elseif N==3 && isrealscalar(x{1}) ...
					&& isrealvector(x{2}) && isrealvector(x{3})
				x{2}=[x{1}-sum(x{2}),x{2}(:).'];
				x{3}=[x{1}-sum(x{3}),x{3}(:).'];
				obj.HeightStart=[x(2),x(3)];
			else
				error('BiErfFit:setHeightStart:InvalidInput',...
					['Input to set the HeightStart as a cell must contain ',...
						'either one real scalar, two real vecotrs, or one scalar ',...
						'and two real vectors.']);
			end
		else
			error('BiErfFit:setHeightStart:InvalidInput',...
				['Input to set the HeightStart must be either a real scalar or ',...
					'a cell containg real numbers.']);
		end
	end

	function obj=set.HeightLow(obj,x)
		if isempty(x)
			obj.HeightLow={};
		elseif isrealscalar(x)
			obj.HeightLow={x,x};
		elseif iscell(x)
			N=numel(x);
			if N==1 && isrealscalar(x{1})
				obj.HeightLow=[x(1),x(1)];
			elseif N==2 && isrealvector(x{1}) && isrealvector(x{2})
				obj.HeightLow={x{1}(:).',x{2}(:).'};
			elseif N==3 && isrealscalar(x{1}) ...
					&& isrealvector(x{2}) && isrealvector(x{3})
				x{2}=[x{1}-sum(x{2}),x{2}(:).'];
				x{3}=[x{1}-sum(x{3}),x{3}(:).'];
				obj.HeightLow=[x(2),x(3)];
			else
				error('BiErfFit:setHeightLow:InvalidInput',...
					['Input to set the HeightLow as a cell must contain ',...
						'either one real scalar, two real vecotrs, or one scalar ',...
						'and two real vectors.']);
			end
		else
			error('BiErfFit:setHeightLow:InvalidInput',...
				['Input to set the HeightLow must be either a real scalar or ',...
					'a cell containg real numbers.']);
		end
	end

	function obj=set.HeightUp(obj,x)
		if isempty(x)
			obj.HeightUp={};
		elseif isrealscalar(x)
			obj.HeightUp={x,x};
		elseif iscell(x)
			N=numel(x);
			if N==1 && isrealscalar(x{1})
				obj.HeightUp=[x(1),x(1)];
			elseif N==2 && isrealvector(x{1}) && isrealvector(x{2})
				obj.HeightUp={x{1}(:).',x{2}(:).'};
			elseif N==3 && isrealscalar(x{1}) ...
					&& isrealvector(x{2}) && isrealvector(x{3})
				x{2}=[x{1}-sum(x{2}),x{2}(:).'];
				x{3}=[x{1}-sum(x{3}),x{3}(:).'];
				obj.HeightUp=[x(2),x(3)];
			else
				error('BiErfFit:setHeightUp:InvalidInput',...
					['Input to set the HeightUp as a cell must contain ',...
						'either one real scalar, two real vecotrs, or one scalar ',...
						'and two real vectors.']);
			end
		else
			error('BiErfFit:setHeightUp:InvalidInput',...
				['Input to set the HeightUp must be either a real scalar or ',...
					'a cell containg real numbers.']);
		end
	end

	function obj=set.WidthStart(obj,x)
		if isempty(x)
			obj.WidthStart={};
		elseif isrealvector(x) && numel(x)==2
				obj.WidthStart={x(1),x(2)};
		elseif iscell(x) && numel(x)==2 ...
				&& isrealvector(x{1}) && isrealvector(x{2})
			obj.WidthStart={x{1}(:).',x{2}(:).'};
		else
			error('BiErfFit:setWidthStart:InvalidInput',...
				['Input to set the WidthStart must be either a two-element ',...
					'vector or a cell containing two real vectors.']);
		end
	end

	function obj=set.WidthLow(obj,x)
		if isempty(x)
			obj.WidthLow={};
		elseif isrealvector(x) && numel(x)==2
				obj.WidthLow={x(1),x(2)};
		elseif iscell(x) && numel(x)==2 ...
				&& isrealvector(x{1}) && isrealvector(x{2})
			obj.WidthLow={x{1}(:).',x{2}(:).'};
		else
			error('BiErfFit:setWidthLow:InvalidInput',...
				['Input to set the WidthLow must be either a two-element ',...
					'vector or a cell containing two real vectors.']);
		end
	end

	function obj=set.WidthUp(obj,x)
		if isempty(x)
			obj.WidthUp={};
		elseif isrealvector(x) && numel(x)==2
				obj.WidthUp={x(1),x(2)};
		elseif iscell(x) && numel(x)==2 ...
				&& isrealvector(x{1}) && isrealvector(x{2})
			obj.WidthUp={x{1}(:).',x{2}(:).'};
		else
			error('BiErfFit:setWidthUp:InvalidInput',...
				['Input to set the WidthUp must be either a two-element ',...
					'vector or a cell containing two real vectors.']);
		end
	end

	function obj=set.BaselineStart(obj,x)
		if isempty(x)
			obj.BaselineStart=[];
		elseif isrealscalar(x)
			obj.BaselineStart=x;
		else
			error('BiErfFit:setBaselineStart:InvalidInput',...
				'Input to set the BaselineStart must be a real scalar.');
		end
	end

	function obj=set.BaselineLow(obj,x)
		if isempty(x)
			obj.BaselineLow=[];
		elseif isrealscalar(x)
			obj.BaselineLow=x;
		else
			error('BiErfFit:setBaselineLow:InvalidInput',...
				'Input to set the BaselineLow must be a real scalar.');
		end
	end

	function obj=set.BaselineUp(obj,x)
		if isempty(x)
			obj.BaselineUp=[];
		elseif isrealscalar(x)
			obj.BaselineUp=x;
		else
			error('BiErfFit:setBaselineUp:InvalidInput',...
				'Input to set the BaselineUp must be a real scalar.');
		end
	end

	function obj=set.Robust(obj,x)
		ME=MException('BiErfFit:setRobust:InvalidInput',...
			'Input to set the Robust must be either: %s.',...
			strjoin(obj.RobustList,', '));
		if isstringscalar(x)
			tf=strcmpi(x,obj.RobustList);
			if any(tf)
				obj.Robust=obj.RobustList{tf};
			else
				throw(ME);
			end
		else
			throw(ME);
		end
	end

	function obj=set.Algorithm(obj,x)
		ME=MException('BiErfFit:setAlgorithm:InvalidInput',...
			'Input to set the Algorithm must be either: %s.',...
			strjoin(obj.AlgorithmList,', '));
		if isstringscalar(x)
			tf=strcmpi(x,obj.AlgorithmList);
			if any(tf)
				obj.Algorithm=obj.AlgorithmList{tf};
			else
				throw(ME);
			end
		else
			throw(ME);
		end
	end

	function obj=set.MovMeanWidth(obj,x)
		ME=MException('BiErfFit:setMovMeanWidth:InvalidInput',...
			['Input to set the MovMeanWidth must be a positive integer scalar ',...
				'or a real scalar between [0,1].']);
		if isrealscalar(x) && x>=0
			if isintegerscalar(x) || x<=1
				obj.MovMeanWidth=x;
			else
				throw(ME);
			end
		else
			throw(ME);
		end
	end

	function obj=set.DiffMaxChange(obj,x)
		if isrealscalar(x) && x>0
			obj.DiffMaxChange=x;
		else
			error('BiErfFit:setDiffMaxChange:InvalidInput',...
				'Input to set the DiffMaxChange must be a positive real scalar.');
		end
	end

	function obj=set.DiffMinChange(obj,x)
		if isrealscalar(x) && x>0
			obj.DiffMinChange=x;
		else
			error('BiErfFit:setDiffMinChange:InvalidInput',...
				'Input to set the DiffMinChange must be a positive real scalar.');
		end
	end

	function obj=set.MaxFunEvals(obj,x)
		if isintegerscalar(x) && x>0
			obj.MaxFunEvals=x;
		else
			error('BiErfFit:setMaxFunEvals:InvalidInput',...
				'Input to set the MaxFunEvals must be a positive integer scalar.');
		end
	end

	function obj=set.MaxIters(obj,x)
		if isintegerscalar(x) && x>0
			obj.MaxIters=x;
		else
			error('BiErfFit:setMaxIters:InvalidInput',...
				'Input to set the MaxIters must be a positive integer scalar.');
		end
	end

	function obj=set.TolFun(obj,x)
		if isrealscalar(x) && x>0
			obj.TolFun=x;
		else
			error('BiErfFit:setTolFun:InvalidInput',...
				'Input to set the TolFun must be a positive real scalar.');
		end
	end

	function obj=set.TolX(obj,x)
		if isrealscalar(x) && x>0
			obj.TolX=x;
		else
			error('BiErfFit:setTolX:InvalidInput',...
				'Input to set the TolX must be a positive real scalar.');
		end
	end

	function obj=set.Center(obj,x)
		if isempty(x)
			obj.Center={};
		elseif iscell(x) && numel(x)==2 ...
				&& isrealmatrix(x{1}) && isrealmatrix(x{2})...
				&& size(x{1},1)==3 && size(x{2},1)==3
			obj.Center=x;
		else
			error('BiErfFit:setCenter:InvalidInput',...
				'Input to set the Center must be a real matrix with 3 rows.');
		end
	end

	function obj=set.Height(obj,x)
		if isempty(x)
			obj.Height={};
		elseif iscell(x) && numel(x)==2 ...
				&& isrealmatrix(x{1}) && isrealmatrix(x{2})...
				&& size(x{1},1)==3 && size(x{2},1)==3
			obj.Height=x;
		else
			error('BiErfFit:setHeight:InvalidInput',...
				'Input to set the Height must be a real matrix with 3 rows.');
		end
	end

	function obj=set.Width(obj,x)
		if isempty(x)
			obj.Width={};
		elseif iscell(x) && numel(x)==2 ...
				&& isrealmatrix(x{1}) && isrealmatrix(x{2})...
				&& size(x{1},1)==3 && size(x{2},1)==3
			obj.Width=x;
		else
			error('BiErfFit:setWidth:InvalidInput',...
				'Input to set the Width must be a real matrix with 3 rows.');
		end
	end

	function obj=set.Baseline(obj,x)
		if isempty(x)
			obj.Baseline=[];
		elseif isrealvector(x) && numel(x)==3
			obj.Baseline=x(:);
		else
			error('BiErfFit:setBaseline:InvalidInput',...
				'Input to set the Baseline must be a real vector of length 3.');
		end
	end

	function obj=set.RelStDev(obj,x)
		if isrealscalar(x) && x>=0
			obj.RelStDev=x;
		else
			error('BiErfFit:setRelStDev:InvalidInput',...
				'Input to set the RelStDev must be a positive real scalar.');
		end
	end

	function obj=set.CoeffDeterm(obj,x)
		if isrealscalar(x)
			obj.CoeffDeterm=x;
		else
			error('BiErfFit:setCoeffDeterm:InvalidInput',...
				'Input to set the CoeffDeterm must be a real scalar.');
		end
	end

	function obj=set.AdjCoeffDeterm(obj,x)
		if isrealscalar(x)
			obj.AdjCoeffDeterm=x;
		else
			error('BiErfFit:setAdjCoeffDeterm:InvalidInput',...
				'Input to set the AdjCoeffDeterm must be a real scalar.');
		end
	end

	function obj=set.NumFunEvals(obj,x)
		if isintegerscalar(x) && x>=0
			obj.NumFunEvals=x;
		else
			error('BiErfFit:setNumFunEvals:InvalidInput',...
				'Input to set the NumFunEvals must be a positive integer scalar.');
		end
	end

	function obj=set.NumIters(obj,x)
		if isintegerscalar(x) && x>=0
			obj.NumIters=x;
		else
			error('BiErfFit:setNumIters:InvalidInput',...
				'Input to set the NumIters must be a positive integer scalar.');
		end
	end

	function obj=set.ExitFlag(obj,x)
		if isempty(x)
			obj.ExitFlag=[];
		elseif isintegerscalar(x)
			obj.ExitFlag=x;
		else
			error('BiErfFit:setExitFlag:InvalidInput',...
				'Input to set the ExitFlag must be an integer scalar.');
		end
	end

	% display the object properties
	disp(obj)

	% construct a fit model from the fit results
	[yModel,yErf,yBaseline]=model(obj,varargin)
end

methods (Access=protected)
	% fit erfs
	obj=fit(obj)
end

methods (Static=true)
	% error function
	y=fnerf(x,varargin)
	
	% Gaussian function
	y=fngaussian(x,c,a,w)
end

end
