function C = combination(n,k)
% function C = combination(n,k)
% Calculates combinates for vector input n and k.  Both inputs need to be
% the same size, or one can be a scalar.
%
% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2010, Eric Tytell <tytell at jhu dot edu>

if (numel(n) == 1),
    sz = size(k);
    n = repmat(n,sz);
elseif (numel(k) == 1),
    sz = size(n);
    k = repmat(k,sz);
elseif ((ndims(n) == ndims(k)) && all(size(n) == size(k))),
    sz = size(n);
else
    error('n and k need to be the same size or scalar.');
end;

N = prod(sz);

n = n(:);
k = k(:);

C = zeros(N,1);
for i = 1:N,
    C(i) = prod((n(i)-k(i)+1):n(i))/factorial(k(i));
end;

C = reshape(C,sz);

