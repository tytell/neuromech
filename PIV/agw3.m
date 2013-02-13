function [up,vp,wp] = agw3(x,y,z, u,v,w, xp,yp,zp, H)

ispivpos = 0;

if (nargin < 10),
  H = [];
  if (nargin < 9),
    xp = [];
    yp = [];
    zp = [];
  end;
end;

if (isempty(xp)),
  xp = x;
  yp = y;
  zp = z;
  ispivpos = 1;
end;

dx = diff(x,[],2);
dy = diff(y,[],1);
dz = diff(z,[],3);

x(isnan(u)) = NaN;

if (all(dx(isfinite(dx)) > 0) & all(dy(isfinite(dy)) > 0)),
  if (any(z ~= repmat(z(1,1,:),[size(z,1) size(z,2)]))),
    error('z value within pages must be the same');
  end;

  dx = nanmean(dx(:));
  dy = nanmean(dy(:));
  dz = nanmean(dz(dz > 0));

  delta = min([dx dy dz]);

  if (isempty(H)),
    H = 1.24*delta;
  end;
  
  % calculate the number of points that actually influence a point in z
  % the effective radius of the window (out to the exp(-3) level)
  rwindz = sqrt(3)*H;
  
  [zp0,q,zpnum] = unique(squeeze(zp(1,1,:)));
  z0 = squeeze(z(1,1,:));

  timedWaitBar(0, 'Smoothing...');

  for i = 1:length(zp0),
    pg1 = find((z0 >= zp0(i)-rwindz) & (z0 <= zp0(i)+rwindz));
    pg2 = find(zpnum == i);

    x1 = flatten(x(:,:,pg1));
    y1 = flatten(y(:,:,pg1));
    z1 = flatten(z(:,:,pg1));
    u1 = flatten(u(:,:,pg1));
    v1 = flatten(v(:,:,pg1));
    w1 = flatten(w(:,:,pg1));
    len1 = length(x1);

    x2 = flatten(xp(:,:,pg2));
    y2 = flatten(yp(:,:,pg2));
    z2 = flatten(zp(:,:,pg2));
    len2 = length(x2);

    x1 = repmat(x1,[1 len2]);
    y1 = repmat(y1,[1 len2]);
    z1 = repmat(z1,[1 len2]);
    u1 = repmat(u1,[1 len2]);
    v1 = repmat(v1,[1 len2]);
    w1 = repmat(w1,[1 len2]);

    x2 = repmat(x2',[len1 1]);
    y2 = repmat(y2',[len1 1]);
    z2 = repmat(z2',[len1 1]);

    dist2 = (x1-x2).^2 + (y1-y2).^2 + (z1-z2).^2;
    mindist = min(dist2);
    k = find(mindist > delta^2);
    dist2(:,k) = NaN;

    alpha = exp(-dist2./H^2);

    aa = nansum(alpha);

    up1 = nansum(alpha.*u1)./aa;
    vp1 = nansum(alpha.*v1)./aa;
    wp1 = nansum(alpha.*w1)./aa;

    up(:,:,pg2) = reshape(up1,[size(xp,1) size(xp,2) length(pg2)]);
    vp(:,:,pg2) = reshape(vp1,[size(xp,1) size(xp,2) length(pg2)]);
    wp(:,:,pg2) = reshape(wp1,[size(xp,1) size(xp,2) length(pg2)]);

    timedWaitBar(i/length(zp0));
  end;
  timedWaitBar(1);
end;


% delta = min([dx dy dz]);
% 
% if (isempty(H)),
%   H = 1.24*delta;
% end;
% 
% % calculate the number of points that actually influence any new point
% % the effective radius of the window (out to the exp(-3) level)
% rwindx = ceil(sqrt(3)*H/dx);
% rwindy = ceil(sqrt(3)*H/dy);
% rwindz = ceil(sqrt(3)*H/dz);
% 
% [jwind,iwind,kwind] = meshgrid(-rwindx:rwindx,-rwindy:rwindy,-rwindz:rwindz);
% q = find(sqrt(iwind.^2/rwindy^2 + jwind.^2/rwindx^2 + kwind.^2/rwindz^2) <= 1);		% ellipsoidal
% iwind = iwind(q);
% jwind = jwind(q);
% kwind = kwind(q);
% 
% xpind = round((xp-min(x(:)))/dx) + 1;
% ypind = round((yp-min(y(:)))/dy) + 1;
% zpind = round((zp-min(z(:)))/dz) + 1;
% xpind(xpind < 1) = 1;			xpind(xpind > size(x,2)) = size(x,2);
% ypind(ypind < 1) = 1;			ypind(ypind > size(x,1)) = size(x,1);
% zpind(zpind < 1) = 1;			zpind(zpind > size(x,3)) = size(x,3);
% 
% npt = prod(size(xpind));
% nwind = prod(size(iwind));
% 
% if (size(xpind,3) == 1),
%   xpind1 = xpind;
%   ypind1 = ypind;
% end;
% 
% for fr = 1:length(t),
%   if (size(xpind1,3) > 1),
%     ypind1 = ypind(:,:,fr);
%     xpind1 = xpind(:,:,fr);
%     xp1 = xp(:,:,fr);
%     yp1 = yp(:,:,fr);
%   end;
%   zpind1 = repmat(zpind(fr),size(xpind1));
%   
%   iall = repmat(iwind(:),[1 npt]) + repmat(ypind1(:)',[nwind 1]);
%   jall = repmat(jwind(:),[1 npt]) + repmat(xpind1(:)',[nwind 1]);
%   kall = repmat(kwind(:),[1 npt]) + repmat(zpind1(:)',[nwind 1]);
%   
%   iall(iall < 1) = NaN;		iall(iall > size(x,1)) = NaN;
%   jall(jall < 1) = NaN;		jall(jall > size(x,2)) = NaN;
%   kall(kall < 1) = NaN;		kall(kall > length(t)) = NaN;
% 
%   q = find(isnan(iall) | isnan(jall) | isnan(kall));
%   
% % same as sub2ind, but without all the slow overrun checks
%   ind2 = iall + (jall-1)*size(x,1);
%   ind2(q) = 1;
%   indrel = iall + ((jall-1) + (kall-fr)*size(x,2))*size(x,1);
%   indrel(q) = 1;
%   
%   z1 = repmat(z(fr
%   if (ispivpos),
%     if (size(x,3) > 1),
%       x1 = x(:,:,fr + (-rwindz:rwindz));
%       y1 = x(:,:,fr + (-rwindz:rwindz));
%       
%       if (size(x,3) > 1),
%         ind3 = iall + ((jall-1) + (kall-1)*size(x,2))*size(x,1);
%         ind3(q) = 1;
%         
%         distx = x(ind3) - repmat(xp1(:)',[nwind 1]);
%         disty = y(ind3) - repmat(yp1(:)',[nwind 1]);
%         distz = z(ind3) - repmat(zp1(:)',[nwind 1]);
%       end;
%       
%       dist2 = distx.^2 + disty.^2 + distz.^2;
%       dist(q) = NaN;
%       
%       alpha = exp(-dist2/H^2);
%       up(:,:,fr) = reshape(nansum(alpha.*u(indall))./nansum(alpha), [size(xp,1) size(xp,2)]);
%       vp(:,:,fr) = reshape(nansum(alpha.*v(indall))./nansum(alpha), [size(xp,1) size(xp,2)]);;
%     end;



	