function m = histmode(x,dim,nbins)
% m = nanmode(x,[dim],[nbins])
% Most common value ignoring NaNs.  For non-repeating values, the "mode" is
% defined as the mean value within the histogram bin with the most values.
% Defaults to 10 bins spanning the range of each column of x.  The number
% of bins can be specified in nbins.  If nbins is Inf or 'exact', then the
% result is truly the most common value, meaning that that individual value
% is repeated the most times.
%
% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2010, Eric Tytell <tytell at jhu dot edu>

if (nargin < 3),
    nbins = 10;
    if (nargin < 2),
        dim = [];
    end;
end;

if (isempty(x)),
    m = NaN;
end;

sz = size(x);
if ((nargin == 1) || isempty(dim)),
    dim = 1;
    if ((ndims(x) == 2) && (size(x,1) == 1)),
        dim = 2;
        x = permute(x,[2 1]);
    end;
elseif (dim ~= 1),
    x = permute(x,[dim 1:dim-1 dim+1:ndims(x)]);
end;

x = reshape(x,[sz(dim) prod(sz([1:dim-1 dim+1:end]))]);

if (isnumeric(nbins) && isfinite(nbins)),
    for i = 1:size(x,2),
        edges = linspace(min(x(:,i)),max(x(:,i)),nbins+1);
        edges(end) = Inf;
        [n,bins] = histc(x(:,i),edges);
        
        [q,ind] = max(n);
        m(i) = mean(x(bins == ind));
    end;
else
    for i = 1:size(x,2),
        m(i) = mode(x(isfinite(x(:,i)),i));
    end;
end;
    
m = reshape(m,[1 sz([1:dim-1 dim+1:end])]);
if (dim ~= 1),
    m = permute(m,[2:dim 1 dim+1:ndims(x)]);
end;
