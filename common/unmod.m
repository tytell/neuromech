function f = unmod(f,m,delta,dim)
% function y = unmod(f,m,delta)
% delta is optional.  Default is m/10

% deal with optional arguments
if (nargin < 4)
	dim = [];
	if (nargin < 3),
		delta = [];
		if (nargin < 2),
			m = 2*pi;
		end;       
	end;
end;

if (isempty(dim))
	[f, dim] = shiftdim(f);
	dim = dim+1;
else
	f = shiftdim(f,dim-1);
end;

if (isempty(delta))
	delta = m/10;
end;

sz = size(f);
rest = prod(sz(2:end));

% step through any dimensions of f beyond the first
for i = 1:rest,
	% only operate on the finite elements
	q = find(isfinite(f(:,i)));
	
	% take their differences
	d = diff(f(q,i));

	% eliminate differences that are close to m
	k = find((abs(m - abs(d)) < delta) & (d > 0));
	d(k) = d(k) - m;
	k = find((abs(m - abs(d)) < delta) & (d < 0));
	d(k) = d(k) + m;

	% use a cumulative sum to get f back
	f(q(2:end),i) = d;
	f(q,i) = cumsum(f(q,i));
end;

f = shiftdim(f,ndims(f)-dim+1);
