function d = deriv(x,y,dim,ord,iseven)
% function d = deriv(x,y,dim,ord,iseven)
% calculates derivatives based on dy(i)/dx = (y(i+1)-y(i-1))/(x(i+1)-x(i-1))
% approximates the end points by a single point difference, rather than the
% central difference
% ord is the order of the derivative (first derivative, second derivative,
% etc)
% iseven specifies whether the spacing is even.  If so, then x can be empty
% or equal to a scalar that specifies the dx value.
%
% Mercurial revision hash: $Revision: 18f43cd9074e $ $Date: 2010/08/10 21:11:58 $
% Copyright (c) 2010, Eric Tytell

if (nargin < 4),
    ord = 1;
    if (nargin < 3)
        dim = [];
        if (nargin < 2)
            error('You must supply at least 2 arguments.');
        end;
    end;
end;

if (nargin < 5),
    if (isempty(x)),
        dx = 1;
        iseven = true;
    elseif (numel(x) == 1)
        iseven = true;
        dx = x;
    else
        iseven = false;
    end;
end;

%for an unspecified dimension, we use the first one -- unless y is a row
%vector, then we operate on the vector
if (isempty(dim)),
    if ((ndims(y) == 2) && (size(y,1) == 1)),
        dim = 2;
    else
        dim = 1;
    end;
end;

%shift the dimension of interest so that it's the first one
szy = size(y);
pmt = [dim 1:dim-1 dim+1:ndims(y)];
y = permute(y,pmt);

if (iseven && (numel(x) > 1)),
    x = permute(x,pmt);
    dx = x(2,:) - x(1,:);
    dx = reshape(dx,[1 prod(szy(pmt(2:end)))]);
elseif (~iseven),
    if ((ndims(x) == 2) && any(size(x) == 1) && (length(x) == size(y,1))),
        %deal with a vector x
        x = repmat(makecol(x),[1 szy(pmt(2:end))]);
    else
        %otherwise permute x like y
        x = permute(x,pmt);
    end;

    if ((ndims(x) ~= ndims(y)) || any(size(x) ~= size(y))),
        error('x and y must be the same size matrices, or x must be a vector');
    end;
end;

np = size(y,1);

%flatten dimensions that aren't interesting
y = reshape(y,[szy(dim) prod(szy(pmt(2:end)))]);
if (~iseven),
    x = reshape(x,[szy(dim) prod(szy(pmt(2:end)))]);
end;

%and do the derivative
d = y;
for i = 1:ord,
    if (iseven),
        d(2:np-1,:) = (d(3:end,:)-d(1:end-2,:));
        d([np 1],:) = (y([end 2],:)-y([end-1 1],:));
    else
        d(2:np-1,:) = (d(3:end,:)-d(1:end-2,:))./(x(3:end,:)-x(1:end-2,:));
        d([np 1],:) = (y([end 2],:)-y([end-1 1],:))./(x([end 2],:)-x([end-1 1],:));
    end;
end;

if (iseven),
    if (numel(dx) == 1),
        d = d / dx^ord;
    else
        dx = repmat(dx,[szy(dim) 1]);
        d = d ./ (dx.^ord);
    end;
end;

%put d make in the permuted shape
d = reshape(d,szy(pmt));
%and undo the permutation
d = permute(d,[2:dim 1 dim+1:length(szy)]);

	