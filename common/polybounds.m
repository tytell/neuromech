function [y,dy] = polybounds(x,b,S,alpha,obs,simul)
% function [y,dy] = polybounds(x,b,S,alpha,obs,simul)
% Calculates (1-alpha)% confidence intervals for regression lines in
% various ways.  Can calculate CIs for prediction or on the line
% itself, and can also calculate them simultaneously or independently
% for the slope and intercept.
%
% Parameters:
%	x - X values at which to return the CI
%	b - Coefficients of the line, returned by polyfit
%	S - Statistics structure, returned by polyfit
%	alpha - Probabality value
%	obs - If 0, calculates CIs for prediction.  If 1, calculates CIs
%		on the line itself
%	simul - If 1, calculates CIs for slope and intercept simultaneously,
%		Otherwise does them separately.
%
% This procedure taken almost verbatim from polytool.m in the Matlab
% Statistics toolbox, version 3.0.

df = S.df;
normr = S.normr;

tcrit = tinv(1 - alpha/2,df);
fcrit = sqrt(length(b) .* finv(1-alpha, length(b), df));

[y,dy] = polyval(b,x,S);

if (obs),
   e = (dy * sqrt(df) / normr).^2;
   e = sqrt(max(0, e-1));
   dy = normr/sqrt(df)*e;
end;

if (simul)
	crit = fcrit;
else
	crit = tcrit;
end;

dy = dy*crit;
