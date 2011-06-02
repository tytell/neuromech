function f = fillmat(a,dim,rev)
%FILLMAT  Fills in NaN values with the previous non-NaN value
%   f = fillmat(a,dim,rev)
%   Fills in NaN values along dimension dim of a with the previous non-NaN
%   value.  If the optional parameter rev is true, starts at the end.
%
% Examples:
%   >> a = [1 2 3 NaN NaN 4 NaN 5];
%   >> fillmat(a,2)
%   ans =
%      1     2     3     3     3     4     4     5
%   >> fillmat(a,2,true)
%   ans =
%      1     2     3     4     4     4     5     5

% Mercurial revision hash: $Revision$ $Date$
% See https://bitbucket.org/tytell/matlab_variables/overview
% Copyright (c) 2011, Eric Tytell <tytell at jhu dot edu>

if (nargin < 3),
    rev = false;
    if (nargin < 2),
        dim = [];
    end;
end;

if (isempty(dim)),
    if ((ndims(a) == 2) && (size(a,1) == 1)),
        dim = 2;
    else
        dim = 1;
    end;
end;

sz = size(a);
ord = [dim 1:dim-1 dim+1:length(sz)];
sz2 = sz(ord);
a = permute(a,ord);
a = reshape(a,[sz2(1) prod(sz2(2:end))]);

if (rev),
    i1 = size(a,1);
    i2 = 1;
    d = -1;
else
    i1 = 1;
    i2 = size(a,1);
    d = 1;
end;

f = zeros(size(a));
f1 = a(1,:);
for i = i1:d:i2,
    good = ~isnan(a(i,:));
    f1(good) = a(i,good);
    f(i,:) = f1;
end;

f = reshape(f,sz2);
f = permute(f,[2:dim 1 dim+1:length(sz)]);
