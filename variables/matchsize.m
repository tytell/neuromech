function x = matchsize(x,szmat)
%MATCHSIZE   Tries to expand a variable to match the size of another
%   x = matchsize(x,y)
% or
%   x = matchsize(x,sz)
% 
%   If x has singleton dimensions where y has multiples, then it uses
%   repmat to increase the size of x along those dimensions.  The second
%   parameter can either be the matrix itself that we're trying to match,
%   or a row vector containing the size of that matrix (ie, size(y)).
%
% Example:
%   If size(y) = [25 7 3] and size(x) = [1 7 3], then
%   a = matchsize(x,y) is equivalent to 
%   a = repmat(x,[25 1 1])
%
% Note: Doesn't actually call REPMAT, which gives it a slight speed boost,
% but it's still a bit slow if you're calling it repeatedly.
%
% See also REPMAT

% Mercurial revision hash: $Revision$ $Date$
% See https://bitbucket.org/tytell/matlab_variables/overview
% Copyright (c) 2011, Eric Tytell <tytell at jhu dot edu>

%if the matrix is a row vector of all integers, assume it represents
%the sizes we want, not a matrix to match
if ((ndims(szmat) == 2) && (size(szmat,1) == 1) && all(floor(szmat) == szmat)),
    sz = szmat;
else
    sz = size(szmat);
end;

if (ndims(x) > length(sz)),
    error('x must be smaller than the matrix to match');
end;
     
%get the size of the first matrix
sz1 = ones(size(sz));
sz1(1:ndims(x)) = size(x);

dim = find((sz1 == 1) & (sz ~= 1));
a = prod(sz(dim));
rest = 1:length(sz1);
rest = rest(~ismember(rest,dim));
if (any(sz(rest) ~= sz1(rest))),
    error('Could not match sizes');
end;
b = prod(sz(rest));

x = permute(x,[dim rest]);
x = reshape(x,[1 b]);
x = x(ones(a,1),:);
x = reshape(x,[sz(dim) sz(rest)]);
x = ipermute(x,[dim rest]);
