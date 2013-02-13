function v = makerow(v)
%MAKEROW   Converts a column vector into a row, if necessary. 
% row = makerow(v)
% If v is a column vector, transpose it to make it a row.  Otherwise leave
% it alone.
%
% See also MAKECOL

% Mercurial revision hash: $Revision$ $Date$
% See https://bitbucket.org/tytell/matlab_variables/overview
% Copyright (c) 2011, Eric Tytell <tytell at jhu dot edu>

if ((ndims(v) == 2) && (size(v,2) == 1)),
    v = v';
end;
