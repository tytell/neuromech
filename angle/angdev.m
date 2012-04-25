function s = angdev(ang,dim)
% ANGDEV - Angular standard deviation
%
%    s = angdev(ang,dim)
%
% Calculates sample angular standard deviation (also returned by ANGMEAN)
% From Fisher 1993, page 32
%
% SEE ALSO
%   angmean

% Mercurial revision hash: $Revision$ $Date$
% See http://rcn.ccs.tulane.edu/index.php5/Tytell_Matlab
% Copyright (c) 2012, Eric Tytell <tytell at jhu dot edu>

if (nargin == 1),
  dim = 1;
end;

angx = nanmean(cos(ang),dim);
angy = nanmean(sin(ang),dim);

r = sqrt(angx.^2 + angy.^2);
lr = log(r);
lr(lr > 0) = NaN;

s = sqrt(-2*lr);
