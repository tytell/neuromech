function sem = nansem(x,dim)
% function sem = nansem(x,dim)
% Like nanstd, but calculates the standard error of the mean
% (stdev / sqrt(n)).
%
% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2010, Eric Tytell <tytell at jhu dot edu>

if (nargin == 1)
    dim = [];
end;

sem = nanstd(x,[],dim) ./ sum(isfinite(x),dim);
