function v = makecol(v)
%MAKECOL   Converts a row vector into a column, if necessary. 
% v = makecol(v)
% If v is a row vector, transpose it to make it a column.  Otherwise leave
% it alone.
%
% See also MAKECOL

% Mercurial revision hash: $Revision$ $Date$
% See https://bitbucket.org/tytell/matlab_variables/overview
% Copyright (c) 2011, Eric Tytell <tytell at jhu dot edu>

if ((ndims(v) == 2) && (size(v,1) == 1)),
    v = v';
end;
