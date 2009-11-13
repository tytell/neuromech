function F = colfilt1(a,w,step, fcn, dim)
% function F = colfilt1(a,w,step, fcn, dim)

if (nargin < 5),
    dim = [];
    if (nargin < 4),
        fcn = @nanmean;
        if (nargin < 3)
            step = 1;
        end;
    end;
end;

if (ischar(step)),
    switch lower(step),
      case 'sliding',
        step = 1;
      case 'distinct',
        step = w;
      otherwise,
        error('Unrecognized step option');
    end;
end;
       
c = vec2col(a,w,step, dim);

F = feval(fcn, c);

if ((ndims(F) > ndims(a)) && (size(F,1) == 1)),
    sz = size(F);
    F = reshape(F,sz(2:end));
elseif ((ndims(a) == 2) && ~isempty(dim) && (dim == 1)),
    F = F';
elseif ((ndims(a) == 2) && (size(a,2) == 1)),
    F = F';
end;

