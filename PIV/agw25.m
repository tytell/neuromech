function [up,vp,wp,H] = agw(x,y,u,v,w,xp,yp,H)
% function [up,vp,wp,H] = agw(x,y,u,v,w,xp,yp,[H])

szu = size(u);
szx = size(x);

if (any(size(x) ~= size(y))),
    error('x and y matrices must be the same size.');
end;
if (any(size(u) ~= size(v))),
    error('u and v matrices must be the same size.');
end;
if (any(szu(1:2) ~= szx(1:2))),
    error('The size of x,y,u, and v must match in the first two dimensions.');
end;
if ((ndims(x) > 2) & ((ndims(x) ~= ndims(u)) | any(szx(3:end) ~= szu(3:end)))),
    error('If x and y are 3D, they must be the same size as u and v.');
end;

if (any(size(xp) ~= size(yp))),
    error('xp and yp must be the same size.');
end;
if (size(xp,3) > size(u,3)),
    error('xp cannot have more frames than u.');
end;

if (nargin < 8),
    H = [];
end;

if (size(xp,3) > 1),
    nfr = size(xp,3);
else
    nfr = size(u,3);
end;

% check to find out if we're approximately plaid
dx = diff(x,[],2);
dy = diff(y,[],1);

% calculate average distances between points
dx = mean(dx(isfinite(dx)));
dy = mean(dy(isfinite(dy)));
dd = sqrt(dx^2 + dy^2);
delta = (2*dx + 2*dy + 4*dd)/8;

nx = size(u,2);
ny = size(u,1);

if (isempty(H)),
    H = 1.24*delta;		% from Agui & Jimenez (1987) JFM 185
end;

% calculate the number of points that actually influence any new point
% the effective radius of the window (out to the exp(-3) level)
rwind = ceil(sqrt(3)*H/delta);

% pad the x,y,u, and v matrices with rwind NaNs
x = [repmat(NaN,[rwind size(x,2)+rwind*2 szx(3:end)]); ...
     repmat(NaN,[size(x,1) rwind szx(3:end)]) x repmat(NaN,[size(x,1) rwind szx(3:end)]); ...
     repmat(NaN,[rwind size(x,2)+rwind*2 szx(3:end)])];
y = [repmat(NaN,[rwind size(y,2)+rwind*2  szx(3:end)]); ...
     repmat(NaN,[size(y,1) rwind szx(3:end)]) y repmat(NaN,[size(y,1) rwind szx(3:end)]); ...
     repmat(NaN,[rwind size(y,2)+rwind*2 szx(3:end)])];
u = [repmat(NaN,[rwind size(u,2)+rwind*2 szu(3:end)]); ...
     repmat(NaN,[size(u,1) rwind szu(3:end)]) u repmat(NaN,[size(u,1) rwind szu(3:end)]); ...
     repmat(NaN,[rwind size(u,2)+rwind*2 szu(3:end)])];
v = [repmat(NaN,[rwind size(v,2)+rwind*2 szu(3:end)]); ...
     repmat(NaN,[size(v,1) rwind szu(3:end)]) v repmat(NaN,[size(v,1) rwind szu(3:end)]); ...
     repmat(NaN,[rwind size(v,2)+rwind*2 szu(3:end)])];
w = [repmat(NaN,[rwind size(w,2)+rwind*2 szu(3:end)]); ...
     repmat(NaN,[size(w,1) rwind szu(3:end)]) w repmat(NaN,[size(w,1) rwind szu(3:end)]); ...
     repmat(NaN,[rwind size(w,2)+rwind*2 szu(3:end)])];

% create the indices of a the points in a circular region rwind in diameter
[ixwind,iywind] = meshgrid(-rwind:rwind);			% square region
q = find(sqrt(ixwind.^2 + iywind.^2) <= rwind);		% circular
ixwind = ixwind(q);
iywind = iywind(q);
n = length(q);

if (nfr > 1),
    timedWaitBar(0,'Smoothing...');
end;
for i = 1:nfr,
    if (size(x,3) > 1),
        x1 = x(:,:,i);
        y1 = y(:,:,i);
    else,
        x1 = x;
        y1 = y;
    end;

    if (size(xp,3) > 1),
        xp1 = xp(:,:,i);
        yp1 = yp(:,:,i);
    else
        xp1 = xp;
        yp1 = yp;
    end;
    
    u1 = u(:,:,i);
    v1 = v(:,:,i);
    w1 = w(:,:,i);

% find the indices in the old field that are closest to each point in the
% new field

    xpindf = ((xp1(:)-min(x1(:)))/dx)' + 1;
    ypindf = ((yp1(:)-min(y1(:)))/dy)' + 1;
    k1 = find((xpindf >= 1) & (xpindf <= nx) & ...
              (ypindf >= 1) & (ypindf <= ny));

    xpindf = xpindf+rwind;
    ypindf = ypindf+rwind;

    closeind(1,:) = floor(ypindf(k1)) + (floor(xpindf(k1))-1)*size(u1,1);
    closeind(2,:) = floor(ypindf(k1)) + (ceil(xpindf(k1))-1)*size(u1,1);
    closeind(3,:) = ceil(ypindf(k1)) + (floor(xpindf(k1))-1)*size(u1,1);
    closeind(4,:) = ceil(ypindf(k1)) + (ceil(xpindf(k1))-1)*size(u1,1);

    k2 = find(any(isfinite(u1(closeind))) & ...
              any(isfinite(v1(closeind))));
    nonan = k1(k2);

    xpind = round(xpindf(nonan));
    ypind = round(ypindf(nonan));
    xp1 = xp1(nonan);
    yp1 = yp1(nonan);

    N = length(nonan);

% when we add the indices of the points closest to the new points to the
% indices of a generic window, we get the indices of the old points
% within rwind points of of each new point
    ix = repmat(xpind,[n 1]) + repmat(ixwind,[1 N]);
    iy = repmat(ypind,[n 1]) + repmat(iywind,[1 N]);

    ind = iy + (ix-1)*size(u1,1);
    
% now calculate the gaussian
    alpha = exp(-((x1(ind) - repmat(xp1(:)',[n 1])).^2 + ...
                  (y1(ind) - repmat(yp1(:)',[n 1])).^2)./H^2);
    uu = alpha.*u1(ind);
    vv = alpha.*v1(ind);
    ww = alpha.*w1(ind);

% and use it to make the new u and v matrices
    up1 = repmat(NaN, [size(xp,1) size(xp,2)]);
    up1(nonan) = nansum(uu)./nansum(alpha);
    vp1 = repmat(NaN, [size(xp,1) size(xp,2)]);
    vp1(nonan) = nansum(vv)./nansum(alpha);
    wp1 = repmat(NaN, [size(xp,1) size(xp,2)]);
    wp1(nonan) = nansum(ww)./nansum(alpha);
    
    up(:,:,i) = up1;
    vp(:,:,i) = vp1;
    wp(:,:,i) = wp1;
    
    if ((nfr > 1) & ~timedWaitBar(i/nfr)),
        break;
    end;
end;

