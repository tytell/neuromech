function d = deriv(x,y,dim)
% function d = deriv(x,y,dim)
% calculates derivatives based on dy(i)/dx = (y(i+1)-y(i-1))/(x(i+1)-x(i-1))
% approximates the end points by a single point difference, rather than the
% central difference

if (nargin < 3)
	dim = [];
	if (nargin < 2)
		error('You must supply at least 2 arguments.');
	end;
end;

if (isempty(dim))
	[y,dim] = shiftdim(y);
	x = shiftdim(x);
	dim = dim+1;
else
	y = shiftdim(y,dim-1);
	x = shiftdim(x,dim-1);
end;

if (size(x,1) ~= size(y,1))
	error('The x and y matrices must have the same number of rows');
end;

if (size(x,2) == 1),
	x = repmat(x,[1 size(y,2)]);
end;
np = size(y,1);

d(2:np-1,:) = (y(3:end,:)-y(1:end-2,:))./(x(3:end,:)-x(1:end-2,:));
d([np 1],:) = (y([end 2],:)-y([end-1 1],:))./(x([end 2],:)-x([end-1 1],:));

d = shiftdim(d,ndims(y)-dim+1);

		