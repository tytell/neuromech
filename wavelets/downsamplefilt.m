function [d,edge] = downsamplefilt(x,n,dim,op)
% function [d,edge] = downsamplefilt(x,n,dim,op)

if (nargin == 4),
    isdim = true;
elseif (nargin == 3),
    if (ischar(dim) || isa(dim,'function_handle')),
        op = dim;
        dim = 1;
        isdim = false;
    else
        isdim = true;
        op = @max;
    end;
elseif (nargin == 2),
    dim = 1;
    op = @max;
    isdim = false;
end;

if (~isdim && (ndims(x) == 1) && (size(x,1) == 1)),
    dim = 2;
end;

if (n - floor(n) ~= 0),
    warning('Downsample factor is not an integer.  Rounding...');
    n = round(n);
end;

sz = size(x);
p = [dim 1:dim-1 dim+1:ndims(x)];
rest = prod(sz(p(2:end)));

x = permute(x,p);
x = reshape(x,[sz(dim) rest]);

len = floor(sz(dim)/n);
edge = floor((sz(dim) - len*n)/2)+1;

x = x(edge+(0:n*len-1),:);
x = reshape(x,[n len rest]);
d = feval(op,x);

d = shiftdim(d,1);
d = reshape(d, [len sz(p(2:end))]);
d = permute(d, [2:dim 1 dim+1:length(sz)]);

edge = edge + ceil(n/2);


function x = middle(x)

m = ceil(size(x,1)/2);
x = x(m,:,:);




