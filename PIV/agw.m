function [up,vp,DU,H] = agw(x,y,u,v,varargin)
% function [up,vp,DU,H] = agw(x,y,u,v,xp,yp,[H])
% Adaptive gaussian window smoothing and interpolation.  Any new point is a sum of
% points nearby, weighted by a Gaussian of window size H.
%
% i.e. up(i) = sum(exp(-((x(j) - xp(i))^2 + (y(j) - yp(i))^2)/H^2 * u(j)) /
% sum(exp(-((x(j) - xp(i))^2 + (y(j) - yp(i))^2)/H^2))
%
% x,y,u, and v are the existing velocity field.  They either have to be
% approximately plaid (i.e. x(i,j) < x(i,j+1), but x(i,j) does not have to
% equal x(i+1,j)), or just vectors of arbitrary positions and velocities.
% For approximately plaid matrices, agw only sums nearby vectors when
% calculating the new value.  For arbitrary position velocity fields, it
% currently calculates a complete distance matrix, from each point in the
% old field to each point in the new one, and uses that to calculate the
% Gaussian.  In the future, if there are many vectors in the arbitrary
% position field, it should iterate through them so that it doesn't have to
% build an enormous distance matrix.  xp and yp are the new positions to for
% which calculate u and v values.  They can be in any form, as long as
% they're the same size.  H is the optional Gaussian window size.  According
% to Agui and Jimenez 87 the optimal value for H is 1.24 delta, where delta
% is the average spacing between vectors in the original field.  If H is not
% passed, agw calculates this value and uses it.
% 
% The approximately plaid specification is required to allow the appropriate
% PIV smoothing operation, which is
%               [us,vs] = agw(x+u/2, y+v/2, u,v, x,y); 
% because the best position for each displacement calculated by PIV is
% actually at the center of the displacement.
%
% If the output is plaid, we calculate a new derivative structure based on
% the smoothed data.

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
if ((ndims(x) > 2) && ((ndims(x) ~= ndims(u)) || any(szx(3:end) ~= szu(3:end)))),
    error('If x and y are 3D, they must be the same size as u and v.');
end;

opt = first(cellfun(@ischar,varargin));
if (isempty(opt)),
    opt = length(varargin)+1;
end;

H = [];
if (opt > 7),
    [xp,yp,H,xr,yr,ur,vr] = deal(varargin{1:7});
    isextra = true;
elseif (opt > 3),
    [xp,yp,H] = deal(varargin{1:3});
    isextra = false;
elseif (opt > 2),
    [xp,yp] = deal(varargin{1:2});
    isextra = false;
    H = [];
else
    xp = x;
    yp = y;

    x = repmat(x,[1 1 size(u,3)]) + u/2;
    y = repmat(y,[1 1 size(v,3)]) + v/2;
    szx = szu;
    
    warning('Assuming that positions and velocities are both is pixels. Is this correct?');
end;

fillnans = false;
isextra = false;

while (opt <= length(varargin)),
    switch lower(varargin{opt}),
        case 'extra',
            [xr,yr,ur,vr] = deal(varargin{opt+(1:4)});
            isextra = true;
            opt = opt + 5;
           
        case 'fillnans',
            fillnans = true;
            opt = opt + 1;
            
        otherwise,
            error('Unrecognized option %s',varargin{opt});
    end;
end;

if (any(size(xp) ~= size(yp))),
    error('xp and yp must be the same size.');
end;
if (size(xp,3) > size(u,3)),
    error('xp cannot have more frames than u.');
end;

if (size(xp,3) > 1),
    nfr = size(xp,3);
else
    nfr = size(u,3);
end;

up = zeros(size(xp,1),size(xp,2),nfr);
vp = zeros(size(xp,1),size(xp,2),nfr);

% check to find out if we're approximately plaid
dx = diff(x,[],2);
dy = diff(y,[],1);

if all(dy(isfinite(dy)) < 0)
    isflipped = true;
    
    x = flipud(x);
    y = flipud(y);
    u = flipud(u);
    v = flipud(v);
else
    isflipped = false;
end

dy = diff(y,[],1);

if (all(size(x)>1) && all(dx(isfinite(dx)) > 0) && all(dy(isfinite(dy)) > 0))
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
    xpad = NaN([szx(1)+2*rwind szx(2)+2*rwind szx(3:end)]);
    ypad = NaN([szx(1)+2*rwind szx(2)+2*rwind szx(3:end)]);
    upad = NaN([szu(1)+2*rwind szu(2)+2*rwind szu(3:end)]);
    vpad = NaN([szu(1)+2*rwind szu(2)+2*rwind szu(3:end)]);
    
    i = rwind + (1:szx(1));
    j = rwind + (1:szx(2));
    xpad(i,j,:) = x;
    ypad(i,j,:) = y;
    upad(i,j,:) = u;
    vpad(i,j,:) = v;

    % create the indices of a the points in a circular region rwind in diameter
    [ixwind,iywind] = meshgrid(-rwind:rwind);			% square region
    iscirc = sqrt(ixwind.^2 + iywind.^2) <= rwind;		% circular
    ixwind = ixwind(iscirc);
    iywind = iywind(iscirc);
    n = sum(iscirc(:));
    
    if (nfr > 1),
        timedWaitBar(0,'Smoothing...');
    end;
    for i = 1:nfr,
        if (size(xpad,3) > 1),
            x1 = xpad(:,:,i);
            y1 = ypad(:,:,i);
        else
            x1 = xpad;
            y1 = ypad;
        end;

        if (size(xp,3) > 1),
            xp1 = xp(:,:,i);
            yp1 = yp(:,:,i);
        else
            xp1 = xp;
            yp1 = yp;
        end;
        
        u1 = upad(:,:,i);
        v1 = vpad(:,:,i);
        
        % find the indices in the old field that are closest to each point in the
        % new field
        
        %get the first finite index, which will include an offset due to
        %rwind, plus if there are any NaN elements around the edges of x1
        %all indices are relative to this one
        xind0 = first(any(isfinite(x1)));
        x0 = first(x1(:,xind0),isfinite(x1(:,xind0)));
        yind0 = first(any(isfinite(y1),2));
        y0 = first(y1(yind0,:),isfinite(y1(yind0,:)));
        %fractional subscripts
        xpindf = ((xp1(:)-x0)/dx)' + xind0;
        ypindf = ((yp1(:)-y0)/dy)' + yind0;
        
        %check to make sure they're in bounds 
        inbounds = (xpindf >= rwind+1) & (xpindf <= size(x1,2)-rwind) & ...
            (ypindf >= rwind+1) & (ypindf <= size(x1,1)-rwind);

        %there are four elements close to two fractional subscripts,
        %depending on whether you round each of them up or down
        %if all of the four elements are NaN, we don't want to estimate a
        %value for that index
        %so find the four elements
        closeind = zeros(4,sum(inbounds));
        closeind(1,:) = floor(ypindf(inbounds)) + (floor(xpindf(inbounds))-1)*size(u1,1);
        closeind(2,:) = floor(ypindf(inbounds)) + (ceil(xpindf(inbounds))-1)*size(u1,1);
        closeind(3,:) = ceil(ypindf(inbounds)) + (floor(xpindf(inbounds))-1)*size(u1,1);
        closeind(4,:) = ceil(ypindf(inbounds)) + (ceil(xpindf(inbounds))-1)*size(u1,1);

        %check to see if any of them are finite
        nonan = false(size(xpindf));
        if (fillnans),
            nonan(inbounds) = true;
        else
            nonan(inbounds) = any(isfinite(u1(closeind))) & any(isfinite(v1(closeind)));
        end;

        %and only keep those indices
        xpind = round(xpindf(nonan));
        ypind = round(ypindf(nonan));
        xp1 = xp1(nonan);
        yp1 = yp1(nonan);

        N = sum(nonan);

        % when we add the subscripts of the points closest to the new points to the
        % subscripts of a generic window, we get the indices of the old points
        % within rwind points of of each new point
        ix = repmat(xpind,[n 1]) + repmat(ixwind,[1 N]);
        iy = repmat(ypind,[n 1]) + repmat(iywind,[1 N]);

        %turn subscripts into indices
        ind = iy + (ix-1)*size(u1,1);
	
        % now calculate the gaussian
        alpha = exp(-((x1(ind) - repmat(xp1(:)',[n 1])).^2 + ...
                      (y1(ind) - repmat(yp1(:)',[n 1])).^2)./H^2);
        uu = alpha.*u1(ind);
        vv = alpha.*v1(ind);
        
        if (isextra),
            xr1 = xr(:,:,i);
            yr1 = yr(:,:,i);
            ur1 = ur(:,:,i);
            vr1 = vr(:,:,i);
            nr = length(xr(:));
            
            dist = sqrt( (repmat(xp1(:)',[nr 1]) - repmat(xr1(:),[1 N])).^2 + ...
                         (repmat(yp1(:)',[nr 1]) - repmat(yr1(:),[1 N])).^2);
            
            alphar = exp(-dist.^2./H^2);
            uur = alphar.*repmat(ur1(:),[1 N]);
            vvr = alphar.*repmat(vr1(:),[1 N]);
            
            alpha = [alpha; alphar];
            uu = [uu; uur];
            vv = [vv; vvr];
        end;

        
        % and use it to make the new u and v matrices
        up1 = NaN([size(xp,1) size(xp,2)]);
        up1(nonan) = nansum(uu)./nansum(alpha);
        vp1 = NaN([size(xp,1) size(xp,2)]);
        vp1(nonan) = nansum(vv)./nansum(alpha);
        
        if ~fillnans && (size(xp,1) == size(x,1)) && ...
                (size(xp,2) == size(x,2))
            good = ~isnan(u(:,:,i));
            up1(~good) = NaN;
            vp1(~good) = NaN;
        end
        
        up(:,:,i) = up1;
        vp(:,:,i) = vp1;
        
        if ((nfr > 1) && ~timedWaitBar(i/nfr)),
            break;
        end;
    end;
else
    warning('Assuming randomly distributed vectors. Is this what you want?');
    
    if (size(u,3) > 1),
        error('Cannot handle 3D non-plaid matrices');
    end;
    
    k = find(isfinite(x) & isfinite(y) & isfinite(u) & isfinite(v));

    x = shiftdim(x(k));
    y = shiftdim(y(k));
    u = shiftdim(u(k));
    v = shiftdim(v(k));

    n = length(x);

    k = find(isfinite(xp) & isfinite(yp));
    np = length(k);
    
    % Assume even, random spacing of points
    delta = sqrt(range(x).*range(y)/n);
    
    if (nargin == 6),
        H = 1.24*delta;
    end;

    % matrix of the distances between each old point (rows) and each new point
    % (columns)
    dist = sqrt( (repmat(xp(k)',[n 1]) - repmat(x,[1 np])).^2 + ...
                 (repmat(yp(k)',[n 1]) - repmat(y,[1 np])).^2);
    
    % Gaussian
    alpha = exp(-dist.^2/H.^2);

    up = NaN(size(xp));
    vp = NaN(size(yp));
    
    % new u and v values
    up(k) = sum(alpha .* repmat(u(:),[1 np]))./sum(alpha);
    vp(k) = sum(alpha .* repmat(v(:),[1 np]))./sum(alpha);
end;

if isflipped
    xp = flipud(xp);
    yp = flipud(yp);
    up = flipud(up);
    vp = flipud(vp);
end

DU.dudx = NaN(size(up));
DU.dudy = NaN(size(up));
DU.dvdx = NaN(size(up));
DU.dvdy = NaN(size(up));

% if the output is plaid, we can recalculate the derivative structure
if ((size(up,1) > 2) && (size(up,2) > 2)),
    if (size(xp,3) > 1),
        DU.dudx(:,2:end-1,:) = (up(:,3:end,:) - up(:,1:end-2,:)) ./ ...
            (xp(:,3:end,:) - xp(:,1:end-2,:));
        DU.dudy(2:end-1,:,:) = (up(3:end,:,:) - up(1:end-2,:,:)) ./ ...
            (yp(3:end,:,:) - yp(1:end-2,:,:));
        DU.dvdx(:,2:end-1,:) = (vp(:,3:end,:) - vp(:,1:end-2,:)) ./ ...
            (xp(:,3:end,:) - xp(:,1:end-2,:));
        DU.dvdy(2:end-1,:,:) = (vp(3:end,:,:) - vp(1:end-2,:,:)) ./ ...
            (yp(3:end,:,:) - yp(1:end-2,:,:));
    else
        DU.dudx(:,2:end-1,:) = (up(:,3:end,:) - up(:,1:end-2,:))./ ...
            repmat(xp(:,3:end) - xp(:,1:end-2), [1 1 size(up,3)]);
        DU.dudy(2:end-1,:,:) = (up(3:end,:,:) - up(1:end-2,:,:))./ ...
            repmat(yp(3:end,:,:) - yp(1:end-2,:,:), [1 1 size(up,3)]);
        DU.dvdx(:,2:end-1,:) = (vp(:,3:end,:) - vp(:,1:end-2,:))./ ...
            repmat(xp(:,3:end,:) - xp(:,1:end-2,:), [1 1 size(up,3)]);
        DU.dvdy(2:end-1,:,:) = (vp(3:end,:,:) - vp(1:end-2,:,:))./ ...
            repmat(yp(3:end,:,:) - yp(1:end-2,:,:), [1 1 size(up,3)]);
    end;
end;

