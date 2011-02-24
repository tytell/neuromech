function C = vec2col(v,w,step, dim)
% function C = vec2col(v,w,step, dim)
% Splits a vector (or the columns of a matrix) into windows with length w,
% stepping with length step.
%
% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2010, Eric Tytell <tytell at jhu dot edu>

if (nargin < 4),
    dim = [];
    if (nargin < 3),
        step = 1;
    end;
end;

sz = size(v);
if (isempty(dim))
    if ((ndims(v) == 2) && (sz(1) == 1)),
        dim = 2;
    else
        dim = 1;
    end;
end;

%permute v so that the dimension of interest is the first
pmt = [dim 1:dim-1 dim+1:ndims(v)];
v = permute(v,pmt);
sz = sz(pmt);
norig = sz(1);
N = prod(sz(2:end));

%flatten the other dimensions
v = reshape(v, [norig N]);

if (isempty(step)),
    step = 1;
end;

%generate the index into v
%first by stepping the center point
ind = 1:step:norig;
nnew = length(ind);

%offsets for the window
off = (0:w-1)' - floor((w-1)/2);

%create the index into the main dimension
ind = ind(ones(w,1),:) + off(:,ones(1,nnew));
good = (ind >= 1) & (ind <= norig);

%replicate for the other dimensions
ind = ind(:,:,ones(1,N));
off = zeros(1,1,N);
off(1,1,:) = (0:N-1)*norig;
off = off(ones(w,1),ones(1,nnew),:);
ind = ind + off;

good = good(:,:,ones(1,N));

%do the indexing
C = NaN(size(ind));
C(good) = v(ind(good));

%and make C the right shape
C = reshape(C, [w nnew sz(2:end)]);
C = ipermute(C, [1 pmt+1]);

%handle row vectors correctly
if ((dim == 2) && (ndims(C) == 3) && (size(C,2) == 1))
    C = reshape(C, [w nnew]);
end;

