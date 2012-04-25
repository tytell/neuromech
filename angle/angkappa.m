function kml = angkappa(ang)
% ANGKAPPA  Calculates dispersion parameter
%
%   kml = angkappa(ang)
%
% Calculates the von Mises dispersion parameter kappa based on maximum
% likelihood estimators.  From Mardia, 2000

% Mercurial revision hash: $Revision$ $Date$
% See http://rcn.ccs.tulane.edu/index.php5/Tytell_Matlab
% Copyright (c) 2012, Eric Tytell <tytell at jhu dot edu>

[mu,r] = angmean(ang);

kml = NaN([1 size(ang,2)]);

r1 = r < 0.53;
r2 = (r >= 0.53) & (r < 0.85);
r3 = (r >= 0.85) & (r < 1);
r4 = r == 1;

kml(r1) = 2*r(r1) + r(r1).^3 + 5*r(r1).^5/6;
kml(r2) = -0.4 + 1.39*r(r2) + 0.43./(1-r(r2));
kml(r3) = 1./(2*(1-r(r3)) - (1-r(r3)).^2 - (1-r(r3)).^3);
kml(r4) = Inf;


