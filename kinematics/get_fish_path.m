function fishpath = get_fish_path(t,x,y,varargin)

opt.smoothper = 1;      % in sec
opt.smoothmethod = 'loess';
opt.ampthresh = 0.01;   % fraction of body length
opt.npt = 20;
opt.path = 'head';  % or 'COM' 
opt.width = ones(opt.npt,1);

if (nargin >= 7) && isnumeric(varargin{1}) && ...
        all(size(varargin{1}) == size(x))
    [hu,hv,hax,hay] = varargin{1:4};
    p = 5;
else
    p = 1;
end
opt = parsevarargin(opt,varargin(p:end),p+3);

dt = t(2)-t(1);
assert(isfinite(dt));
if isvector(x)
    good = isfinite(x) & isfinite(y);
else
    good = any(isfinite(x) & isfinite(y));
end
smoothfrac = opt.smoothper / (sum(good)*dt);

switch opt.path
    case 'head'
        good = isfinite(x) & isfinite(y);
        pathx = NaN(size(x));
        pathy = NaN(size(y));

        if (smoothfrac < 1)
            pathx(good) = smooth(x(good),smoothfrac,opt.smoothmethod);
            pathy(good) = smooth(y(good),smoothfrac,opt.smoothmethod);
        else
            p = polyfit(t,x(good),2);
            pathx(good) = polyval(t(good),p);
            p = polyfit(t,y(good),2);
            pathy(good) = polyval(t(good),p);
        end

        pathcurve = curvature(pathx,pathy, 'smooth',1, 'splineindiv');

        swimvecx = [diff(pathx) NaN];
        swimvecy = [diff(pathy) NaN];
        mag = sqrt(swimvecx.^2 + swimvecy.^2);
        swimvecx = swimvecx ./ mag;
        swimvecy = swimvecy ./ mag;

        pathang = unwrap(atan2(swimvecy,swimvecx)+pi) - pi;

        speed = hu.*swimvecx + hv.*swimvecy;
        accel = hax.*swimvecx + hay.*swimvecy;

    case {'COM','com'}
        %crude estimate of COM position - weighted average of midline
        %position, weighted by body width
        comx = NaN(1,size(x,2));
        comy = NaN(1,size(x,2));
        comx(good) = nansum(bsxfun(@times,x(:,good), opt.width)) ./ nansum(opt.width);
        comy(good) = nansum(bsxfun(@times,y(:,good), opt.width)) ./ nansum(opt.width);
        
        pathx = comx;
        pathy = comy;
        
        npts = sum(isfinite(x));
        
        offx = 0;
        offy = 0;
        for i = 2:length(comx)
            if good(i) && good(i-1) && (npts(i) ~= npts(i-1))
                offx = -comx(i) + pathx(i-1);
                offy = -comy(i) + pathy(i-1);
            end
            pathx(i) = pathx(i) + offx;
            pathy(i) = pathy(i) + offy;
        end
                
        pathx(good) = smooth(pathx(good),smoothfrac,opt.smoothmethod);
        pathy(good) = smooth(pathy(good),smoothfrac,opt.smoothmethod);

        isstep = [false npts(2:end) ~= npts(1:end-1)];
        pathx(isstep) = NaN;
        pathy(isstep) = NaN;
        
        comu = deriv(t,pathx);
        comv = deriv(t,pathy);
        isstep(1:end-1) = isstep(1:end-1) | isstep(2:end);
        comu(isstep) = NaN;
        comv(isstep) = NaN;
        
        comax = deriv(t,comu);
        comay = deriv(t,comv);
        
        pathcurve = curvature(pathx,pathy, 'smooth',0, 'discrete');
        
        swimvecx = [diff(pathx) NaN];
        swimvecy = [diff(pathy) NaN];
        mag = sqrt(swimvecx.^2 + swimvecy.^2);
        swimvecx = swimvecx ./ mag;
        swimvecy = swimvecy ./ mag;

        pathang = unwrap(atan2(swimvecy,swimvecx)+pi) - pi;

        speed = comu.*swimvecx + comv.*swimvecy;
        accel = comax.*swimvecx + comay.*swimvecy;
        
        fishpath.comx = comx;
        fishpath.comy = comy;
        fishpath.comu = comu;
        fishpath.comv = comv;
        fishpath.comax = comax;
        fishpath.comay = comay;
end

fishpath.pathx = pathx;
fishpath.pathy = pathy;
fishpath.pathcurve = pathcurve;
fishpath.pathang = pathang;
fishpath.swimvecx = swimvecx;
fishpath.swimvecy = swimvecy;
fishpath.speed = speed;
fishpath.accel = accel;


