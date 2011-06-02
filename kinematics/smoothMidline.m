function [mxs,mys,oxs,oys] = smoothMidline(t,mx,my,len,width, method,varargin)
% function [mxs,mys,oxs,oys] = smoothMidline(t,mx,my,len,width, method,varargin)
% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2010, Eric Tytell

nfr = size(mx,2);
npt = size(mx,1);

if (any((sum(isfinite(mx)) > 1) & (sum(isfinite(mx)) < npt))),
    s = repmat(NaN,size(mx));
    for i = 1:size(mx,2),
        k = find(isfinite(mx(:,i)) & isfinite(my(:,i)));
        if (~isempty(k)),
            s(k,i) = [0; cumsum(sqrt(diff(mx(k,i)).^2 + ...
                                     diff(my(k,i)).^2))];
        end;
    end;
else
    s = [zeros(1,nfr); cumsum(sqrt(diff(mx).^2 + diff(my).^2))];
end;

fprintf('Digitized length %.2f+-%.2f = %.1f%% of input (%.2f)\n', ...
        nanmean(s(end,:)), nanstd(s(end,:))/sqrt(sum(isfinite(s(end,:)))), ...
        nanmean(s(end,:))/len*100, len);

if ((abs((nanmean(s(end,:))-len)./len) > 0.05) | ...
		(nanstd(s(end,:))./len  > 0.05)),
    warning('Lengths don''t match up (either digitized length varies or is different from input).');
end;

switch lower(method),
 case 'spaps',
  serr = varargin{1};
  terr = varargin{2};
  if (length(varargin) == 3),
      goodpts = varargin{3};
  else
      goodpts = 1:npt;
  end;

  ds = diff(s);
  ds0 = len/(npt-1);
  if (any((ds - ds0)/ds0 > 0.05)),
      error('Distance along midline increases weirdly.');
  end;

  s0 = nanmedian(s,2);
  s = linspace(0,len,npt);

  k = find(all(isfinite(mx) & isfinite(my)));
  a = k(1):k(end);

  XY = cat(1,shiftdim(mx(goodpts,k),-1),shiftdim(my(goodpts,k),-1));
  sp = spaps({s0(goodpts),t(k)}, XY, ...
             {serr^2*range(s0(goodpts))*length(t(k)) ...
              terr^2*range(t(k))*length(s0(goodpts))});
  XYs = fnval(sp,{s,t(a)});

  mxs = repmat(NaN,size(mx));
  mys = repmat(NaN,size(my));
  mxs(:,a) = squeeze(XYs(1,:,:));
  mys(:,a) = squeeze(XYs(2,:,:));

  dxyds = fnval(fnder(sp,[1 0]),{s,t(a)});


 case 'spapsang',
  serr = varargin{1};
  terr = varargin{2};

  ang = atan2(diff(my),diff(mx));
  
  s0 = nanmedian(s,2);
  s = linspace(0,len,npt);
  seglen = len/(npt-1);

  k = find(all(isfinite(mx) & isfinite(my)));
  a = k(1):k(end);

  % fudge the length slightly
  if (s0(end-1) < s(end-1)),
      s0(end-1) = s(end-1);
  end;
  sp = spaps({s0(1:end-1),t(k)}, ang(:,k), ...
             {serr^2*range(s0(1:end-1))*length(t(k)) ...
              terr^2*range(t(k))*(length(s0)-1)});
  angs = fnval(sp,{s(1:end-1),t(a)});

  dxyds(1,2:npt-1,a) = (cos(angs(1:end-1,:)) + cos(angs(2:end,:)))/2;
  dxyds(1,[1 npt],a) = cos(angs([1 end],:));
  dxyds(2,2:npt-1,a) = (sin(angs(1:end-1,:)) + sin(angs(2:end,:)))/2;
  dxyds(2,[1 npt],a) = sin(angs([1 end],:));
  dxyds(:,:,[1:k(1)-1 k(end)+1:nfr]) = NaN;

  mxs = zeros(size(mx));
  mys = zeros(size(my));
  mxs(1,:) = mx(1,:);
  mys(1,:) = my(1,:);
  mxs(2:end,a) = repmat(mx(1,a),[npt-1 1]) + cumsum(cos(angs)*seglen);
  mys(2:end,a) = repmat(my(1,a),[npt-1 1]) + cumsum(sin(angs)*seglen);
  mxs(:,[1:k(1)-1 k(end)+1:nfr]) = NaN;
  mys(:,[1:k(1)-1 k(end)+1:nfr]) = NaN;

 case 'agw',
  good = varargin{1};
  ratio = varargin{2};

  k = find(any(isfinite(mx) & isfinite(my)));
  a = k(1):k(end);

  tt = repmat(t,[npt 1]);
  dt = mean(diff(t(k)));
  tt = tt/dt;

  ds = nanmean(flatten(diff(s(:,k))));
  s = s/ds/ratio;

  ss = repmat(linspace(0,len,npt)',[1 nfr])/ds/ratio;

  mx(~good) = NaN;
  my(~good) = NaN;

  [mxs,mys,DU] = agw(tt(:,k),s(:,k),mx(:,k),my(:,k), tt(:,a),ss(:,a), ...
                     []);

  dxyds(1,:,:) = DU.dudy;
  dxyds(2,:,:) = DU.dvdy;
end;

width = shiftdim(width);
if (length(width) ~= npt),
    width = spline(linspace(0,1,length(width)),width, ...
                   linspace(0,1,npt))';
end;
oxs = repmat(NaN,[2*npt nfr]);
oys = repmat(NaN,[2*npt nfr]);
w = repmat(width,[1 k(end)-k(1)+1])/2 * len;
oxs(1:npt,a) = mxs(:,a) - w .* squeeze(dxyds(2,:,a));
oxs(2*npt:-1:npt+1,a) = mxs(:,a) + ...
    w .* squeeze(dxyds(2,:,a));
oys(1:npt,a) = mys(:,a) + w .* squeeze(dxyds(1,:,a));
oys(2*npt:-1:npt+1,a) = mys(:,a) - ...
    w .* squeeze(dxyds(1,:,a));

