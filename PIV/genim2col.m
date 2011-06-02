function t = genim2col(A, varargin)
% function 	t = genim2col(A, x,y, nbh)
%	or 		t = genim2col(A, x,y, [m n])
%	or		t = genim2col(A, nbh)
%	or		t = genim2col(A, [m n])
%
% Generalized version of im2col.  Converts the matrix A into a set of regions
% in the neighborhood of the points in x and y.  For each point (x,y), the corresponding
% column of t are the points within a rectangular region in A, m rows and n columns,
% centered on (x,y).  For a more complex region (circular, for example), you can use nbh
% to specify which of the neighboring points are included.  In the rectangular region
% of the same size as nbh, points where nbh are 1 are included in t.  The central point
% of nbh is ceil(size(nbh)/2).
%
% If you don't specify the x and y points, then genim2col operates on each point in A.
%
% Any points that fall beyond the edge of the matrix A are given NaN values in t.

% Check the arguments
nbh = [];
x = [];
y = [];
if (nargin == 2),
	if (prod(size(varargin{1})) == 2),
		sz = varargin{1};
		m = sz(1);
		n = sz(2);
	else
		nbh = varargin{1};
	end;
elseif (nargin == 4),
	x = varargin{1};
	y = varargin{2};
	if (prod(size(varargin{3})) == 2),
		sz = varargin{3};
		m = sz(1);
		n = sz(2);
	else
		nbh = varargin{3};
	end;
else
	error('Wrong number of arguments.');
end;

% If no (x,y) are passed, then we use all points
if (isempty(x)),
	[x,y] = meshgrid(1:size(A,2), 1:size(A,1));
end;

% Make the points row vectors
x = x(:)';
y = y(:)';
npt = length(x);

% sizes
[M,N] = size(A);
if (~isempty(nbh)),
	[m,n] = size(nbh);
end;
m2 = ceil(m/2);
n2 = ceil(n/2);

[nx,ny] = meshgrid((1:n) - n2, (1:m) - m2);
if (~isempty(nbh)),
	nx = nx(nbh == 1);		% Only keep points where nbh is 1
	ny = ny(nbh == 1);
else
	nx = nx(:);
	ny = ny(:);
end;
Nnbh = length(nx);		% number of points in the neighborhood

% create x and y coordinates of neighborhood points
tx = repmat(x,[Nnbh 1]) + repmat(nx,[1 npt]);
ty = repmat(y,[Nnbh 1]) + repmat(ny,[1 npt]);
% flag those that are off the image
offim = (tx < 1) | (tx > N) | (ty < 1) | (ty > N);

% create indices for the x and y coordinates
% NB: we don't use sub2ind here, because it does a lot of checks that we don't
% need, which is slow
ind = (tx-1).*M + ty;
ind(offim) = 1;				% use a dummy index for the off image points

t = A(ind);					% get the points
t(offim) = NaN;

