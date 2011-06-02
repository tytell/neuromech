function e = appxeq(a,b,tol)
%APPXEQ   Test for two values equal to within a percent tolerance
%   e = appxeq(a,b,tol)
%   a and b are approximately equal if their difference is within tol
%   percent of their mean.

% Mercurial revision hash: $Revision$ $Date$
% See https://bitbucket.org/tytell/matlab_variables/overview
% Copyright (c) 2011, Eric Tytell <tytell at jhu dot edu>

if (nargin == 2),
    tol = 0.05;
end;

e = abs((a-b)./(0.5*(abs(a)+abs(b)))) < tol;
