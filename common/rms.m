function r = rms(a, dim)
% r = rms(a, dim)
% Calculates the RMS value across a particular dimension.
%
% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2010, Eric Tytell <tytell at jhu dot edu>

%handle when they don't pass in a dimension
if (nargin == 1),
    %for a vector, operate along dimension 1 or 2, depending on which one
    %has the data
    if (ndims(a) == 2),
        dim = first(size(a) ~= 1);
    else
        dim = 1;
    end;
end;

%rearrange a so that the dimension of interest is first, and everything
%else is flattened
sz = size(a);
nd = ndims(a);
notdim = [1:dim-1 dim+1:nd];
a = permute(a,[dim notdim]);
a = reshape(a,[sz(dim) prod(sz(notdim))]);

%now get rid of NaN values
good = isfinite(a);
a(~good) = 0;

%squared
a2 = a.^2;

%root mean
r = sqrt(sum(a2) ./ sum(good));

%rearrange r so that it matches what a was
r = reshape(r,[1 sz(notdim)]);
r = permute(r,[2:dim 1 dim+1:nd]);
