function disp(obj)
%% Display Object Properties
%
% Copyright: Herianto Lim
% http://heriantolim.com/
% First created: 29/05/2013
% Last modified: 25/10/2016

fprintf('\nFit Options:\n');
fprintf('\t%24s:    %s\n','Method',obj.Method);
fprintf('\t%24s:    %s\n','Robust',obj.Robust);
fprintf('\t%24s:    %s\n','Algorithm',obj.Algorithm);
fprintf('\t%24s:    %.2e\n','DiffMaxChange',obj.DiffMaxChange);
fprintf('\t%24s:    %.2e\n','DiffMinChange',obj.DiffMinChange);
fprintf('\t%24s:    %d\n','MaxFunEvals',obj.MaxFunEvals);
fprintf('\t%24s:    %d\n','MaxIters',obj.MaxIters);
fprintf('\t%24s:    %.2e\n','TolFun',obj.TolFun);
fprintf('\t%24s:    %.2e\n','TolX',obj.TolX);

exitFlag=obj.ExitFlag;
if ~isempty(exitFlag)
	numErfs=obj.NumErfs;
	header1='   Lower Bound   Start Point   Upper Bound';
	header2='      Lower CI   Convergence      Upper CI     RelUncert\n';
	format1='   %11.4e   %11.4e   %11.4e';
	format2='   %11.4e   %11.4e   %11.4e   %11.4e\n';
	
	fprintf('\nCenter:\n\tCurve   Erf');
	fprintf(header1);
	fprintf(header2);
	lower=obj.CenterLow;
	start=obj.CenterStart;
	upper=obj.CenterUp;
	x=obj.Center;
	for i=1:2
		for j=1:numErfs(i)
			fprintf('\t%5d   %3d',i,j);
			fprintf(format1,lower{i}(j),start{i}(j),upper{i}(j));
			fprintf(format2,x{i}(2,j),x{i}(1,j),x{i}(3,j),...
				abs((x{i}(3,j)-x{i}(2,j))/x{i}(1,j))/2);
		end
	end
	
	fprintf('\nHeight:\n\tCurve   Erf');
	fprintf(header1);
	fprintf(header2);
	lower=obj.HeightLow;
	start=obj.HeightStart;
	upper=obj.HeightUp;
	x=obj.Height;
	for i=1:2
		for j=1:numErfs(i)
			fprintf('\t%5d   %3d',i,j);
			fprintf(format1,lower{i}(j),start{i}(j),upper{i}(j));
			fprintf(format2,x{i}(2,j),x{i}(1,j),x{i}(3,j),...
				abs((x{i}(3,j)-x{i}(2,j))/x{i}(1,j))/2);
		end
	end
	
	fprintf('\nWidth:\n\tCurve   Erf');
	fprintf(header1);
	fprintf(header2);
	lower=obj.WidthLow;
	start=obj.WidthStart;
	upper=obj.WidthUp;
	x=obj.Width;
	for i=1:2
		for j=1:numErfs(i)
			fprintf('\t%5d   %3d',i,j);
			fprintf(format1,lower{i}(j),start{i}(j),upper{i}(j));
			fprintf(format2,x{i}(2,j),x{i}(1,j),x{i}(3,j),...
				abs((x{i}(3,j)-x{i}(2,j))/x{i}(1,j))/2);
		end
	end
	
	fprintf('\nBaseline:\n\t           ');
	fprintf(header1);
	fprintf(header2);
	lower=obj.BaselineLow;
	start=obj.BaselineStart;
	upper=obj.BaselineUp;
	x=obj.Baseline;
	for i=1:2
		for j=1:numErfs(i)
			fprintf('\t%5d   %3d',i,j);
			fprintf(format1,lower(j),start(j),upper(j));
			fprintf(format2,x(2,j),x(1,j),x(3,j),abs((x(3,j)-x(2,j))/x(1,j))/2);
		end
	end
	
	fprintf('\nHysteresis Profile:\n\t%52s ','');
	fprintf(header2);
	x=obj.Asymptote;
	fprintf('\t%52s:','Asymptote(:,1)');
	fprintf(format2,x(2,1),x(1,1),x(3,1),abs((x(3,1)-x(2,1))/x(1,1))/2);
	fprintf('\t%52s:','Asymptote(:,2)');
	fprintf(format2,x(2,2),x(1,2),x(3,2),abs((x(3,2)-x(2,2))/x(1,2))/2);
	x=obj.DynamicRange;
	fprintf('\t%52s:','DynamicRange');
	fprintf(format2,x(2,1),x(1,1),x(3,1),abs((x(3,1)-x(2,1))/x(1,1))/2);
	x=obj.HysteresisWidth;
	fprintf('\t%52s:','HysteresisWidth');
	fprintf(format2,x(2,1),x(1,1),x(3,1),abs((x(3,1)-x(2,1))/x(1,1))/2);
	x=obj.TransitionWidth;
	for i=1:2
		fprintf('\t%52s:',sprintf('TransitionWidth{%d}',i));
		fprintf(format2,x{i}(2,1),x{i}(1,1),x{i}(3,1),...
			abs((x{i}(3,1)-x{i}(2,1))/x{i}(1,1))/2);
	end
	
	fprintf('\nGoodness of Fit:\n');
	fprintf('\t%24s:    %.4g\n','RelStDev',obj.RelStDev);
	fprintf('\t%24s:    %.4g\n','CoeffDeterm',obj.CoeffDeterm);
	fprintf('\t%24s:    %.4g\n','AdjCoeffDeterm',obj.AdjCoeffDeterm);

	fprintf('\nIteration Report:\n');
	fprintf('\t%24s:    %d\n','NumFunEvals',obj.NumFunEvals);
	fprintf('\t%24s:    %d\n','NumIters',obj.NumIters);
	if exitFlag>0
		exitStatus='converged';
	elseif exitFlag<0
		exitStatus='diverged';
	elseif obj.NumFunEvals>obj.MaxFunEvals
		exitStatus='MaxFunEvals was exceeded';
	elseif obj.NumIters>obj.MaxIters
		exitStatus='MaxIters was exceeded';
	else
		exitStatus='fail';
	end
	fprintf('\t%24s:    %d (%s)\n','ExitFlag',exitFlag,exitStatus);
end

end