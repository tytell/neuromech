function [xs,ys,us,vs,axs,ays] = interps(x,y, err)
% function [xs,ys,us,vs,axs,ays] = interps(x,y, err)
% Smooth and interpolate the x and y data using a smoothing spline.  Data
% is smoothed so that the mean squared error between the smoothed and
% original data in the 2D Cartesian space is approximately equal to err.
% Also estimates velocity and acceleration.
%
% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2010, Eric Tytell <tytell at jhu dot edu>

if (nargin == 2),
	err = 0.5;
end;

if (any(size(x) ~= size(y))),
	error('x and y must be the same size');
end;

[x,n] = shiftdim(x);
y = shiftdim(y);

fr = 1:length(x);

k = find(isfinite(x) & isfinite(y));

sp = spaps(fr(k), [x(k) y(k)]', err^2*range(fr(k)), 3);
xys = fnval(sp, fr);
xys(:,[1:k(1)-1 k(end)+1:end]) = NaN;

xs = shiftdim(xys(1,:)',-n);
ys = shiftdim(xys(2,:)',-n);

uvs = fnval(fnder(sp), fr);
uvs(:,[1:k(1)-1 k(end)+1:end]) = NaN;

us = shiftdim(uvs(1,:)',-n);
vs = shiftdim(uvs(2,:)',-n);

axys = fnval(fnder(sp,2), fr);
axys(:,[1:k(1)-1 k(end)+1:end]) = NaN;

axs = shiftdim(axys(1,:)',-n);
ays = shiftdim(axys(2,:)',-n);

