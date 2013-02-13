function cpivPlotError(varargin)

if (nargin >= 3),
	x = varargin{1};
	y = varargin{2};
	error = varargin{3};
	p = 4;
else
	error = varargin{1};
	[x,y] = meshgrid(1:size(error,2),1:size(error,1));
	p = 2;
end;

if (nargin == p),
	scale = varargin{p};
else
	scale = 1;
end;

k = find(error > 0);
x = x(k);
y = y(k);
error = error(k);

dx = [-1 1];
dy = [-1; 1];

u = zeros(size(x));
v = zeros(size(y));

u(bitget(error,1) == 1) = dx(1);
u(bitget(error,2) == 1) = dx(2);
v(bitget(error,3) == 1) = dy(1);
v(bitget(error,4) == 1) = dy(2);

rest = bitget(error,5) + bitget(error,6);

quiverc(x,y,u,v,rest,scale);
