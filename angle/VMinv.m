function theta = VMinv(P, mu,kappa)
% function theta = VMinv(P, mu,kappa)
% Inverse of the von Mises CDF.  Returns angles theta for a set of
% probabilities P, given a mean and concentration of mu and kappa.
%
% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2010, Eric Tytell


N = numel(P);

theta0 = linspace(0,2*pi,100);
P0 = VMcdf(theta0, mu,kappa);

theta1 = csapi(P0,theta0, P);

P2 = VMcdf(theta1, mu,kappa);

if (any(abs(P2-P) > 1e-6)),
  warning('VMinv:badfit','Bad fit to VM distribution');
end;

theta = theta1;

