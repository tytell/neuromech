function m = flatten(m,d) 
%FLATTEN   Flattens dimensions of a matrix
%  m = flatten(m,d)
%  Flattens certain dimensions of m.  The dimensions to flatten, d, is
%  optional.  If it's not provided, flatten is equivalent to m(:).  If it is
%  provided, the dimensions specified in d are merged into the first
%  dimension of d.  Inspired by a similar command in Mathematica.
%
% Example:
%   size(m) = [5 8 3 2]
%   size(flatten(m)) = [5*8*3*2 1]
%   size(flatten(m,[3 4])) = [5 8 6]
%   size(flatten(m,[1 4])) = [10 8 3]

% Mercurial revision hash: $Revision$ $Date$
% See https://bitbucket.org/tytell/matlab_variables/overview
% Copyright (c) 2011, Eric Tytell <tytell at jhu dot edu>


if (nargin == 1),
    m = m(:);
else
    nd = ndims(m);
    sz = size(m);
    if (max(d) > nd),
        sz(end+1:max(d)) = 1;
        nd = max(d);
    end;
    nond = setdiff(1:nd,d);

    m = permute(m,[nond d]);
    if (length(d) == nd),
        m = m(:);
    else
        m = reshape(m,[sz(nond) prod(sz(d))]);
        nd = ndims(m);

        pm = 1:nd;
        pm(d(1)+1:end) = pm(d(1):end-1);
        pm(d(1)) = nd;
        m = permute(m,pm);
    end;
end;

