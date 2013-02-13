function F = VMcdf(theta, mu,kappa)
% function F = VMcdf(theta, mu,kappa)
% Calculates the CDF for the Von-Mises distribution with parameters mu and
% kappa, from 0 to theta.
%
% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2010, Eric Tytell

N = numel(theta);
F = zeros(size(theta));
for i = 1:N,
  if (theta(i) ~= 0),
    F(i) = quad(@VMcf, 0, theta(i), [],[], mu,kappa);
  end;
end;

F = F / (2*pi*besseli(0,kappa));
F = reshape(F,size(theta));



function a = VMcf(phi, mu,kappa)

a = exp(kappa*cos(phi - mu));
