function agwArb(x0,y0,u,v, xp,yp, H)

delta = sqrt(range(x0)*range(y0)/prod(size(x0)));

if (nargin == 6),
	H = 1.24*delta;
end;

tri = delaunay(x0,y0);
k0 = dsearch(x0,y0, tri, xp,yp);

