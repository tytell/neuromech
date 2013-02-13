function [sem,conf] = angconf(ang,alpha)
% ANGCONF - Standard error and confidence interval for mean angle
%
%    [sem,conf] = angconf(ang,alpha)
%
%Calculates the parametric standard error and confidence interval for the
%mean angle of ang, assuming a von Mises distribution.
%From Fisher 1993, p. 89

% Mercurial revision hash: $Revision$ $Date$
% See http://rcn.ccs.tulane.edu/index.php5/Tytell_Matlab
% Copyright (c) 2012, Eric Tytell <tytell at jhu dot edu>

if (nargin == 1),
    alpha = 0.05;
end;

[amean,R] = angmean(ang);
kappa = angkappa(ang);
n = sum(isfinite(ang));

%estimate SEM
sem = 1./sqrt(n.*R.*kappa);

%normal distribution critical point
z = norminv(1-0.5*alpha, 0,1);
z = repmat(z,[1 size(ang,2)]);

%confidence interval if the data are too dispersed is just the whole circle
conf(1:2,1:size(ang,2)) = [1;1]*amean + ...
    repmat([-pi/2; pi/2],[1 length(amean)]);

%for better data, we have the formula below
k = find(z.*sem <= 1);
if (~isempty(k)),
    conf(:,k) = [1;1]*amean(k) + [-1;1]*asin(z(k).*sem(k));
end;
