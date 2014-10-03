function varargout = curvature(x,y,varargin)
% function curve = curvature(x,y,dim,method)
% Calculates the curvature of a set of (x,y) curves.
% Method can be
%   'discrete' - Estimates segment angle and takes a numerical derivative
%       of angle with respect to arc length
%   'spline' - Splines the curves using a quintic spline and estimates the curvature using
%       the equation (x'*y'' - y'*x'') / (x'^2 + y'^2)^1.5.  A smoothing parameter (spaps
%       "tol" parameter) can be specified with the 'smooth' option.  The spline is
%       returned as a second output.  The spline is bivariate (ie, in arc length and
%       time).  Assumes that points are evenly spaced along the arc length.
%   'splineindiv' - Splines each frame individually.  Otherwise the same as 'spline', but
%       works when points are missing or aren't evenly distributed along the body length.
%       This method is *much* slower than 'spline'.
%
% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2010, Eric Tytell

opt.dim = [];
opt.method = 'discrete';
opt.smooth = 0;
opt.filterorder = 4;

if ((length(varargin) >= 1) && isnumeric(varargin{1}) && (numel(varargin{1}) == 1)),
    opt.dim = varargin{1};
    i = 2;
else
    i = 1;
end;

opt = parsevarargin(opt,varargin(i:end),'firstoptionnumber',2+i,...
                    'multival',{'method',{'discrete','spline','splineindiv'}});

if ((ndims(x) ~= ndims(y)) || any(size(x) ~= size(y))),
    error('x and y must be the same size');
end;

%check to see if we have a row vector
if (isempty(opt.dim)),
    if (size(x,1) == 1)
        opt.dim = 2;                         
    else
        opt.dim = 1;
    end;
end;

%rearrange the dimensions in x and y so that the dimension we're interested
%in is the first one, and everything else is flattened
sz = size(x);
perm = [opt.dim 1:opt.dim-1 opt.dim+1:length(sz)];
x = permute(x,perm);
x = flatten(x,2:ndims(x));
y = permute(y,perm);
y = flatten(y,2:ndims(y));

switch lower(opt.method),
  case 'discrete',
    %calculate the angle at each node (mean angle of the two segments
    %on either side of the node for interior nodes, angle of the
    %neighboring segment for the two end nodes)
    segang = atan2(diff(y),diff(x));

    %make sure angles are continuous
    segang = unwrap(segang);

    %smooth, if necessary
    if (length(opt.smooth) == 2)
        segang = smooth2d(segang, opt.filterorder,opt.smooth(1),opt.smooth(2));
    elseif ((length(opt.smooth) == 1) && (opt.smooth > 0)),
        [b,a] = butter(opt.filterorder, opt.smooth);
        
        for i = 1:size(segang,2),
            segang1 = filtfilt(b,a, segang(:,i));
            segang(:,i) = segang1;
        end;
    end;
    
    ds = sqrt(diff(x).^2 + diff(y).^2);
    
    %take the derivative
    curve = NaN(size(x));
    curve(2:end-1,:) = diff(segang) ./ ((ds(1:end-1,:)+ds(2:end,:))/2);

    switch nargout,
      case 1,
        varargout = {curve};
        
      case 2,
        s = [zeros(1,size(x,2)); cumsum(ds)];
        varargout = {s,curve};
        
      case 3,
        s = [zeros(1,size(x,2)); cumsum(ds)];
        segang(size(x,1),:) = NaN;      % get size the same as x
        varargout = {s,curve, segang};
    end;
    
  case 'splineindiv',       
    ds = sqrt(diff(x).^2 + diff(y).^2);
    ds(isnan(ds)) = 0;
    s = [zeros(1,size(x,2)); cumsum(ds)];
    s(isnan(x) | isnan(y)) = NaN;
    
    %generate a set of splines for the curves
    %bivariate spline (sp = {x(s),y(s)})
    %quintic along s
    %no smoothing
    XY = cat(1,shiftdim(x,-1),shiftdim(y,-1));
    dxyds = NaN(size(XY));
    dxyds2 = NaN(size(XY));
    
    ind = first(sum(isfinite(s),1) > 3);
    good = all(isfinite(XY(:,:,ind))) & isfinite(s(:,ind)');
    
    if (numel(opt.smooth) == 1),
        opt.smooth = opt.smooth*ones(1,1,size(XY,3));
    elseif (numel(opt.smooth) ~= size(XY,3)),
        error(['Smooth parameter must be scalar or have the same number of frames ' ...
               'as x and y']);
    end;
    
    sp = spaps(s(good,ind),XY(:,good,ind),opt.smooth(1)^2,3);
    sp = repmat(sp,[1 size(x,2)]);
    for i = 1:size(x,2),
        good = all(isfinite(XY(:,:,i))) & isfinite(s(:,i)');
        if (sum(good) >= 3),
            sp(i) = spaps(s(good,i),XY(:,good,i),opt.smooth(i)^2*range(s(good,i)),3);

            %estimate derivatives
            dxyds(:,:,i) = fnval(fnder(sp(i),1),s(:,i));
            dxyds2(:,:,i) = fnval(fnder(sp(i),2),s(:,i));
        end;
    end;

    %curvature
    %for formula, see e.g. Wikipedia
    %http://en.wikipedia.org/wiki/Curvature
    curve = (dxyds(1,:,:).*dxyds2(2,:,:) - dxyds(2,:,:).*dxyds2(1,:,:)) ./ ...
            (dxyds(1,:,:).^2 + dxyds(2,:,:).^2).^1.5;
    
    %the first dimension will have size 1, so get rid of it
    curve = reshape(curve,size(x));
    curve = reshape(curve,sz(perm));
    curve = ipermute(curve,perm);
    
    switch nargout,
      case 1,
        varargout = {curve};
      case 2,
        varargout = {curve,sp};
      case 3,
        varargout = {s,curve,sp};
    end;
        
  case 'spline',
    s = [zeros(1,size(x,2)); cumsum(sqrt(diff(x).^2 + diff(y).^2))];
    fr = 1:size(x,2);
    
    %assume that s is roughly constant (ie, points are evenly spaced and remain so in
    %each frame)
    s = nanmedian(s,2);
    
    %generate a set of splines for the curves
    %bivariate spline (sp = {x(s),y(s)})
    %quintic along s
    %no smoothing
    XY = cat(1,shiftdim(x,-1),shiftdim(y,-1));
    
    good = all(all(isfinite(XY)));

    sp = spaps({s, fr(good)}, XY(:,:,good), {opt.smooth^2*range(s)*sum(good) 0},{3 1});
    
    dxyds = fnval(fnder(sp,[1 0]), {s, fr(good)});
    dxyds2 = fnval(fnder(sp,[2 0]), {s, fr(good)});
    
    %curvature
    %for formula, see e.g. Wikipedia
    %http://en.wikipedia.org/wiki/Curvature
    curve = NaN(1,length(s),length(fr));
    curve(:,:,good) = (dxyds(1,:,:).*dxyds2(2,:,:) - dxyds(2,:,:).*dxyds2(1,:,:)) ./ ...
            (dxyds(1,:,:).^2 + dxyds(2,:,:).^2).^1.5;
    
    %the first dimension will have size 1, so get rid of it
    curve = reshape(curve,size(x));
    
    curve = reshape(curve,sz(perm));
    curve = ipermute(curve,perm);

    switch nargout,
      case 1,
        varargout = {curve};
      case 2,
        varargout = {curve,sp};
      case 3,
        varargout = {s,curve,sp};
    end;
    
  otherwise,
    error('Unrecognized method %s', method);
end;

