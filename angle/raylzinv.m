function P = raylzinv(Z,n)
%RAYLZINV   Computes the inverse of the Rayleigh Z distribution
%
%    P = raylzinv(Z,n)
%
% SEE ALSO
%   ANGMEANTEST

% Mercurial revision hash: $Revision$ $Date$
% See http://rcn.ccs.tulane.edu/index.php5/Tytell_Matlab
% Copyright (c) 2012, Eric Tytell <tytell at jhu dot edu>

if (any(n(:) < 10)),
  warning('Accuracy will be low with low n.');
end;

R = sqrt(Z.*n);
P = exp(sqrt(1+4*n+4*(n.^2-R.^2)) - (1+2*n)); % from Zar, p 617


