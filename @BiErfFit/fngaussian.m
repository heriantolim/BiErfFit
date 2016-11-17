function y=fngaussian(x,c,a,w)
%% Gaussian function with center c, area a, and FWHM w
%
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 04/04/2013
% Last modified: 04/04/2013

y=2/sqrt(pi/log(2))*a/w*exp(-4*log(2)*(x-c).^2/w^2);

end