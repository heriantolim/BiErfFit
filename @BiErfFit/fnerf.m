function y=fnerf(x,varargin)
%% Error Function
%  Note: the erf is a well defined function, there will be no division by zero
%  when width=0, instead the function will become a step function.
%
% Requires package:
%  - Common_v1.0.0+
%
% Tested on:
%  - MATLAB R2015b
%
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 12/06/2013
% Last modified: 25/10/2016

assert(isrealvector(x),...
	'BiErfFit:fnerf:InvalidInput',...
	'Input to the domains must be a vector of real numbers.');
assert(all(cellfun(@isrealscalar,varargin)),...
	'BiErfFit:fnerf:InvalidInput',...
	'Input to the function parameters must be a real scalar.');

N=nargin;
switch N
	case 3
		c=varargin{1};
		h=1;
		w=varargin{2};
		b=0;
	case 4
		c=varargin{1};
		h=varargin{2};
		w=varargin{3};
		b=0;
	case 5
		c=varargin{1};
		h=varargin{2};
		w=varargin{3};
		b=varargin{4};
	otherwise
		error('BiErfFit:fnerf:WrongNargin',...
			'Unexpected number of input arguments.');
end

y=b+h/2*(1+erf(2*sqrt(log(2))*(x-c)/w));

end
