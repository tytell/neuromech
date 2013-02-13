function [a,b] = tformdispinv(tform,x,y,u,v)
% Transforms a displacement vector.  For correct scaling, u and v should be
% in the same units as x and y (i.e., mm displacement for mm coordinates,
% not mm/s)

good = isfinite(u) & isfinite(v);
if ((size(u,3) > 1) && (size(x,3) == 1)),
    x = repmat(x,[1 1 size(u,3)]);
    y = repmat(y,[1 1 size(v,3)]);
end;

xd1 = [x(good) x(good)+u(good)];
yd1 = [y(good) y(good)+v(good)];

[xd2,yd2] = tforminv(tform,xd1,yd1);
a = repmat(NaN,size(u));
b = repmat(NaN,size(v));
a(good) = diff(xd2,[],2);
b(good) = diff(yd2,[],2);
