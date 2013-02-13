function [mx,my, ox,oy] = getMidline3(avi,frames, hx,hy, tx,ty, npts, width,len, varargin)

info = aviinfo(avi);

if (isempty(frames)),
    frames = first(isfinite(hx)):last(isfinite(hx));
end;
nfr = info.NumFrames;

if ((length(hx) < frames(end)) || any(~isfinite(hx(frames))) || ...
        (length(hy) < frames(end)) || any(~isfinite(hy(frames))) || ...
        (length(tx) < frames(end)) || any(~isfinite(tx(frames))) || ...
        (length(ty) < frames(end)) || any(~isfinite(ty(frames)))),
    error('Need head and tail position in each frame');
end;

invert = false;
debug = false;
quiet = false;
eqhist = false;
mx0 = [];
my0 = [];
weight = [4 0.5 2 0.5 2 1 1];
norm = [1 0.08 0.05 6 300 1 1];
cutoff = [0.3 10 10 10 10 2 2];
fps = 50;

npeaks = 2;
curvenear = 1;

i = 1;
while (i <= length(varargin)),
    switch lower(varargin{i}),
        case 'invert',
            if ((i+1 <= length(varargin)) && islogical(varargin{i+1})),
                invert = varargin{i+1};
                i = i+2;
            else
                invert = true;
                i = i+1;
            end;
        case 'eqhist',
            eqhist = true;
            i = i+1;
        case 'startmidline',
            mx0 = varargin{i+1};
            my0 = varargin{i+2};
            i = i+3;
        case 'fps',
            fps = varargin{i+1};
            i = i+2;
        case 'norm',
            norm = varargin{i+1};
            i = i+2;
        case 'cutoff',
            cutoff = varargin{i+1};
            i = i+2;
            
        case 'debug',
            debug = true;
            i = i+1;
        case 'quiet',
            quiet = true;
            i = i+1;
            
        otherwise,
            error('Unrecognized option %s', varargin{i});
    end;
end;

if (~quiet),
    set(gcf,'DoubleBuffer','on');
end;

%width is specified in body lengths, so multiply by length (in pixels) to
%get width in pixels
width = width.*len;
if (length(width) ~= npts),
    width = spline(linspace(0,len,length(width)),width, ...
        linspace(0,len,npts));
end;
width(width < 1) = 1;
width = makecol(width);

maxwidth = max(width);
maxtrans = 2*maxwidth;

%transect coordinates
perps0 = -maxtrans:maxtrans;
perps0 = perps0 + (maxtrans-perps0(end))/2;       %center the transect
perps = repmat(perps0,[npts-2 1]);
nperp = size(perps0,2);

seglen0 = len / (npts-1);
s0 = (0:seglen0:len)';

%generate a first straight midline
if (isempty(mx0)),
    mx1 = linspace(hx(frames(1)),tx(frames(1)),npts)';
    my1 = linspace(hy(frames(1)),ty(frames(1)),npts)';
else
    mx1 = mx0;
    my1 = my0;
end;

norm(4) = norm(4)/fps;
norm(5) = norm(5)/fps^2;
logcutoff = [log(cutoff(1)) -log(cutoff(2:end))];

mx = zeros(npts,length(frames));
my = zeros(npts,length(frames));
ox = zeros(2*npts,length(frames));
oy = zeros(2*npts,length(frames));

prevdx = [];
prevdy = [];
prevds = [];
prevda = [];

der = zeros(npts-2,nperp);

reader = mmreader(avi);

if (~quiet),
    I = im2double(read(reader, frames(1)));
    if (size(I,3) > 1),
        I = rgb2gray(I);
    end;    
    if (invert),
        I = 1-I;
    end;
    
    himage = imshow6(I,'n');
    hmid = addplot(mx1,my1,'r.-');
    hout = addplot([mx1; mx1(end:-1:1)],[my1; my1(end:-1:1)],'g-');
end;

for i = 1:length(frames),
    fr = frames(i);
    set(gcf, 'Name', sprintf('Frame %d', fr));
    
    I = im2double(read(reader, frames(i)));
    if (size(I,3) > 1),
        I = rgb2gray(I);
    end;
    
    if (invert),
        I = 1-I;
    end;
    if (eqhist),
        I = histeq(I);
    end;
    
    %shift the last frame's midline to the current head and tail position
    [mx1,my1] = rotateAndScale(mx1,my1,1,hx(fr),hy(fr),npts,tx(fr),ty(fr));
    s1 = [0; cumsum(sqrt(diff(mx1).^2 + diff(my1).^2))];
    dxds = (mx1(3:end)-mx1(1:end-2))./(s1(3:end) - s1(1:end-2));
    dyds = (my1(3:end)-my1(1:end-2))./(s1(3:end) - s1(1:end-2));
    
    %collect edges
    perpx = repmat(mx1(2:end-1),[1 nperp]) + repmat(dyds,[1 nperp]).*perps;
    perpy = repmat(my1(2:end-1),[1 nperp]) - repmat(dxds,[1 nperp]).*perps;
    
    %intersections between neighboring transects
    dist = zeros(npts-2,2);
    for m = 1:npts-2,
        if (m > 1),
            n = m-1;
            dist(m,1) = (1/dyds(n)*(mx1(n+1)-mx1(m+1)) + ...
                1/dxds(n)*(my1(n+1)-my1(m+1))) / ...
                (dyds(m)/dyds(n) - dxds(m)/dxds(n));
        end;
        if (m < npts-2),
            n = m+1;
            dist(m,2) = (1/dyds(n)*(mx1(n+1)-mx1(m+1)) + ...
                1/dxds(n)*(my1(n+1)-my1(m+1))) / ...
                (dyds(m)/dyds(n) - dxds(m)/dxds(n));
        end;
    end;
    good = (perps.*sign(repmat(dist(:,1),[1 nperp])) <= abs(repmat(dist(:,1),[1 nperp]))) & ...
        (perps.*sign(repmat(dist(:,2),[1 nperp])) <= abs(repmat(dist(:,2),[1 nperp])));
    
    trans = interp2(I,perpx,perpy, '*nearest');
    
    %fourth order numerical derivative: see
    %http://en.wikipedia.org/wiki/Numerical_differentiation
    der(:,3:end-2) = (-trans(:,5:end) + 8*trans(:,4:end-1) - ...
        8*trans(:,2:end-3) + trans(:,1:end-4))/12;
    der(~good) = 0;
    
    %look for left side peaks (negative)
    isleft = false(size(der));
    isleft(:,2:end-1) = (der(:,2:end-1) < 0) & (der(:,2:end-1) < der(:,1:end-2)) & ...
        (der(:,2:end-1) < der(:,3:end));
    %and right side (positive)
    isright = false(size(der));
    isright(:,2:end-1) = (der(:,2:end-1) > 0) & (der(:,2:end-1) > der(:,1:end-2)) & ...
        (der(:,2:end-1) > der(:,3:end));
    
    %now sort them
    lpeaks = der;
    lpeaks(~isleft) = 0;
    [lpeaks,l] = sort(lpeaks,2,'ascend');
    %only keep the largest npeaks
    l = perps0(l(:,1:npeaks));
    lpeaks = lpeaks(:,1:npeaks);
    
    rpeaks = der;
    rpeaks(~isright) = 0;
    [rpeaks,r] = sort(rpeaks,2,'descend');
    r = perps0(r(:,1:npeaks));
    rpeaks = rpeaks(:,1:npeaks);
    
    %then we have npeaks^2 number of corresponding midlines
    nmid = npeaks^2;
    altmx = repmat(mx1,[1 nmid]);
    altmy = repmat(my1,[1 nmid]);
    lx = zeros(npts-2,nmid);
    ly = zeros(npts-2,nmid);
    rx = zeros(npts-2,nmid);
    ry = zeros(npts-2,nmid);
    
    sharp = nans(npts-2,nmid);
    
    a = 1;
    for m = 1:npeaks,
        for n = 1:npeaks,
            %generate the left and right side coordinates from each of the peaks
            good = l(:,m) < r(:,n);
            
            xx = mx1(2:end-1);
            xx(~good) = NaN;
            yy = my1(2:end-1);
            yy(~good) = NaN;
            
            lx(:,a) = xx + dyds .* l(:,m);
            ly(:,a) = yy - dxds .* l(:,m);
            rx(:,a) = xx + dyds .* r(:,n);
            ry(:,a) = yy - dxds .* r(:,n);
            
            altmx(2:end-1,a) = (lx(:,a) + rx(:,a))/2;
            altmy(2:end-1,a) = (ly(:,a) + ry(:,a))/2;

            sharp(good,a) = rpeaks(good,m) - lpeaks(good,n);

            a = a+1;
        end;
    end;
    
    bad = isnan(sharp(:,1));
    if (any(bad)),
        [nextpeak,pk] = max(sharp(bad,:),[],2);
        
        pt = find(bad);
        
        ind = sub2ind(size(lx),pt,pk);
        lx(pt,1) = lx(ind);
        lx(ind) = NaN;
        ly(pt,1) = ly(ind);
        ly(ind) = NaN;
        rx(pt,1) = rx(ind);
        rx(ind) = NaN;
        ry(pt,1) = ry(ind);
        ry(ind) = NaN;
        
        altmx(pt+1,1) = (lx(pt,1) + rx(pt,1))/2;
        altmy(pt+1,1) = (ly(pt,1) + ry(pt,1))/2;
        
        sharp(bad,1) = nextpeak;
        sharp(ind) = NaN;
    end;

    dx = diff(altmx(:,1));
    dy = diff(altmy(:,1));
    ds = sqrt(dx.^2 + dy.^2);
    
    dangp = nans(npts-2,nmid);
    % dot product to give the angle between two segments
    dangp(:,1) = acos((dx(1:end-1).*dx(2:end) + dy(1:end-1).*dy(2:end)) ./ ...
        (ds(1:end-1).*ds(2:end)));
    
    dangp2 = nans(npts-2,nmid);
    dangp2(2:end-1,1) = 2*dangp(2:end-1,1) - dangp(1:end-2,1) - dangp(3:end,1);
    
    dangt = nans(npts-2,nmid);
    if (~isempty(prevdx)),
        da = acos((dx.*prevdx + dy.*prevdy) ./ (ds .* prevds));
        dangt(:,1) = (abs(da(1:end-1)) + abs(da(2:end)))/2;
    end;
    
    dangt2 = nans(npts-2,nmid);
    if (~isempty(prevda)),
        da2 = da - prevda;
        dangt2(:,1) = (abs(da2(1:end-1)) + abs(da2(2:end)))/2;
    end;

    dwidth = abs(sqrt((lx-rx).^2 + (ly-ry).^2) - repmat(width(2:end-1),[1 nmid]));

    seglen = nans(npts-2,nmid);
    ds = sqrt(dx.^2 + dy.^2);
    seglen(:,1) = (ds(1:end-1) + ds(2:end))/2;
    
    %now go through and sub in the secondary peaks one by one and estimate
    %the changes in angle
    for j = 2:npts-1,
        near = j-2:j+2;

        %make sure we don't go over the ends
        good = (near > 0) & (near <= npts);
        near = near(good);

        ctr = find(near == j);
        
        for k = 2:nmid,
            if (isfinite(altmx(j,k))),
                %estimate changes in curvature
                %get the x and y coordinates and estimate derivatives
                xx = altmx(near,1);
                yy = altmy(near,1);

                xx(ctr) = altmx(j,k);
                yy(ctr) = altmy(j,k);

                dx = diff(xx);
                dy = diff(yy);
                ds = sqrt(dx.^2 + dy.^2);
                
                if (~isempty(prevdx)),
                    da = acos((dx.*prevdx(near(1:end-1)) + dy.*prevdy(near(1:end-1))) ./ ...
                        (ds .* prevds(near(1:end-1))));
                    dangt(j-1,k) = abs(da(ctr-1)) + abs(da(ctr));
                end;
                
                if (~isempty(prevda)),
                    da2 = da - prevda(near(1:end-1));
                    dangt2(j-1,k) = abs(da2(ctr-1)) + abs(da(ctr));
                end;
                
                da = acos((dx(1:end-1).*dx(2:end) + dy(1:end-1).*dy(2:end)) ./ ...
                    (ds(1:end-1).*ds(2:end)));
    
                dangp(j-1,k) = da(ctr-1);
                
                if (length(near) == 5),
                    dangp2(j-1,k) = 2*da(ctr-1) - da(ctr-2) - da(ctr);
                end;
                
                %also the change in segment length
                if (ctr > 1),
                    seglen(j-1,k) = (ds(ctr-1)+ds(ctr))/2;
                else
                    seglen(j-1,k) = ds(ctr);
                end;
            end;
        end;
    end;

    %normalize
    sharp = sharp ./ mean(sharp(:,1)) / norm(1);
    dangp = abs(dangp) / norm(2);
    dangp2 = abs(dangp2) / norm(3);
    dangt = abs(dangt) / norm(4);
    dangt2 = abs(dangt2) / norm(5);
    seglen = abs(seglen - seglen0) / seglen0 / norm(6);
    dwidth = dwidth ./ maxwidth / norm(7);
    
    %now maximize sharp while minimizing some function of dpos, dvel, seglen, and 
    %dcurve
    opt = cat(3,log(sharp),-log(dangp),-log(dangp2),-log(dangt),-log(dangt2),...
        -log(seglen),-log(dwidth));
    
    %cutoff points that are beyond the threshold for the different
    %parameters
    bad = opt < matchsize(logcutoff,opt);
    
    logM = nansum(matchsize(weight,opt).*opt, 3);
    logM(all(isnan(opt),3)) = NaN;
    logM(any(bad,3)) = -Inf;
    
    %sort the effects of changing each point individually
    
    %first find the max benefit/cost for each point individually
    [logM0,pk0] = max(logM,[],2);
    ind0 = sub2ind([npts nmid],(2:npts-1)',pk0);
    
    logM1 = logM;
    logM1(sub2ind([npts-2 nmid],(1:npts-2)',pk0)) = -Inf;
    logM1(isnan(logM1)) = -Inf;
    
    %then sort the remaining
    [logMs,ord] = sort(logM1(:), 'descend');
    [pt,pk] = ind2sub([npts-2 nmid],ord);
    
    %start from the midline that results from the sharpest peaks
    testmx = [hx(fr); altmx(ind0); tx(fr)];
    testmy = [hy(fr); altmy(ind0); ty(fr)];
    
    %change the points that have the biggest effects sequentially and step
    %until the value of logM is no longer increasing
    testgoodpk = isfinite(logM0);
    logM0 = sum(logM0(testgoodpk));
    logM1 = logM0;
    prevlogM = -Inf;
    testpk = pk0;
    curpk = testpk;
    j = 1;
    while ((logM1 >= logM0) && (logM1 > prevlogM)),
        %save the current configuration
        mx2 = testmx;
        my2 = testmy;
        curpk = testpk;
        prevlogM = logM1;
        goodpk = testgoodpk;
        
        dx = diff(testmx(:,1));
        dy = diff(testmy(:,1));
        ds = sqrt(dx.^2 + dy.^2);
        
        dp1 = acos((dx(1:end-1).*dx(2:end) + dy(1:end-1).*dy(2:end)) ./ ...
            (ds(1:end-1).*ds(2:end)));
    
        dp2 = nans(npts-2,1);
        dp2(2:end-1) = 2*dp1(2:end-1) - dp1(1:end-2) - dp1(3:end);

        if (~isempty(prevdx)),
            da = acos((dx.*prevdx + dy.*prevdy) ./ ...
                (ds .* prevds));
            dt1 = (abs(da(1:end-1)) + abs(da(2:end)))/2;
        else
            dt1 = nans(size(dp1));
        end;

        if (~isempty(prevda)),
            da2 = da - prevda;
            dt2(:,1) = (abs(da2(1:end-1)) + abs(da2(2:end)))/2;
        else
            dt2 = nans(size(dp1));
        end;

        dp1 = abs(dp1) / norm(2);
        dp2 = abs(dp2) / norm(3);
        dt1 = abs(dt1) / norm(4);
        dt2 = abs(dt2) / norm(5);
        
        seglen = abs((ds(1:end-1)+ds(2:end))/2 - seglen0) / seglen0;
                           
        %and the new value of logM
        ind = sub2ind([npts-2 nmid],(1:npts-2)',testpk);
        opt = cat(3,log(sharp(ind)),-log(dp1),-log(dp2),-log(dt1),-log(dt2),...
            -log(seglen),-log(dwidth(ind)));
        
        bad = opt < matchsize(logcutoff,opt);
        
        logM1 = nansum(matchsize(weight,opt).*opt, 3);
        
        testgoodpk = ~any(bad,3);
        logM1 = sum(logM1(testgoodpk));

        %substitute in the next point and peak (remember they're ordered
        %from biggest effect on logM to smallest)
        testpk(pt(j)) = pk(j);
        ind = sub2ind([npts nmid],(2:npts-1)',testpk);
        
        testmx(2:end-1) = altmx(ind);
        testmy(2:end-1) = altmy(ind);
        j = j+1;
    end;
    
    %interpolate back onto evenly spaced points
    good = isfinite(mx2) & isfinite(my2) & [true; goodpk; true];
    s1 = [0; cumsum(sqrt(diff(mx2(good)).^2 + diff(my2(good)).^2))];
    s0 = linspace(0,s1(end),npts)';
    
    sp = csape(s1',[mx2(good) my2(good)]','variational');
    xy = fnval(sp, s0');
    mx1 = xy(1,:)';
    my1 = xy(2,:)';

    mx(:,fr) = mx1;
    my(:,fr) = my1;

    %estimate the outline
    dxy = fnval(fnder(sp),s0')';
   
    ind = sub2ind(size(lx),(1:npts-2)',curpk);
    width1 = [1; sqrt((lx(ind)-rx(ind)).^2 + (ly(ind)-ry(ind)).^2); 1];
    width1 = interp1(s1,width1(good), s0(2:end-1), 'spline');

    ox([1 npts npts+1 2*npts],fr) = mx2([1 npts npts 1]);
    oy([1 npts npts+1 2*npts],fr) = my2([1 npts npts 1]);
    
    ox(2:npts-1,fr) = mx1(2:npts-1) + dxy(2:end-1,2).*width1/2;
    ox((npts-1:-1:2) + npts,fr) = mx1(2:npts-1) - dxy(2:end-1,2).*width1/2;
    oy(2:npts-1,fr) = my1(2:npts-1) - dxy(2:end-1,1).*width1/2;
    oy((npts-1:-1:2) + npts,fr) = my1(2:npts-1) + dxy(2:end-1,1).*width1/2;
    
    dx = diff(mx1);
    dy = diff(my1);
    ds = sqrt(dx.^2 + dy.^2);

    if (i > 1),
        prevda = acos((dx.*prevdx + dy.*prevdy) ./ ...
            (ds .* prevds));
    end;

    prevdx = dx;
    prevdy = dy;
    prevds = ds;
    
    if (~quiet),
        set(himage,'CData',I);
        set(hmid,'XData',mx2(good),'YData',my2(good));
        set(hout,'XData',ox(:,fr),'YData',oy(:,fr));
        drawnow;
    end;
end;

            
% -------------------------------------------------------
function [x,y] = rotateAndScale(x,y, a,xa,ya, b,xb,yb)

dx1 = x(b)-x(a);
dy1 = y(b)-y(a);
dx2 = xb-xa;
dy2 = yb-ya;
olddist = sqrt(dx1^2 + dy1^2);
newdist = sqrt(dx2^2 + dy2^2);

dotprod = (dx1*dx2 + dy1*dy2) / (olddist*newdist);
sgn = sign(dx1*dy2 - dx2*dy1);

if (dotprod > 1),                       % rounding error problem
    dotprod = 1;
end;

ang = sgn*acos(dotprod);

x2 = x;
y2 = y;
x2(a:b) = ((x(a:b) - x(a))*cos(ang) - (y(a:b) - y(a))*sin(ang)) / ...
      olddist * newdist + xa;
y2(a:b) = ((x(a:b) - x(a))*sin(ang) + (y(a:b) - y(a))*cos(ang)) / ...
      olddist * newdist + ya;

x = x2;
y = y2;

if ((abs(x2(b)-xb) > 1) || (abs(y2(b)-yb) > 1)),
    warning('problems...');
end;

function [ds,dev] = deviation(x,y)

good = isfinite(x) & isfinite(y);

xx = x(good);
yy = y(good);

dx1 = xx(2:end) - xx(1:end-1);
dy1 = yy(2:end) - yy(1:end-1);
ds1 = sqrt(dx1.^2 + dy1.^2);

dx1 = dx1(1:end-1);
dy1 = dy1(1:end-1);

dx2 = xx(3:end) - xx(1:end-2);
dy2 = yy(3:end) - yy(1:end-2);

dev1 = (dx1.*dy2 - dy1.*dx2) ./ sqrt(dx2.^2 + dy2.^2);

dev = nans(size(x));
dev(good) = [NaN; dev1; NaN];

ds = zeros(size(x));
ds(good) = [ds1; 0];


function [ds,cv] = curve1(x,y)

good = isfinite(x) & isfinite(y);

ds1 = sqrt(diff(x(good)).^2 + diff(y(good)).^2);
ds = zeros(size(x));
ds(good) = [ds1; 0];

s = [0; cumsum(ds1)];

dxds = deriv(s,x(good));
dyds = deriv(s,y(good));
dxds2 = deriv(s,dxds);
dyds2 = deriv(s,dyds);

cv = zeros(size(x));
cv(good) = (dxds.*dyds2 - dyds.*dxds2) ./ ...
    ((dxds.^2 + dyds.^2).^1.5);

gapstart = find(good(1:end-1) & ~good(2:end));
gapend = find(~good(1:end-1) & good(2:end));

if (~isempty(gapstart) && ~isempty(gapend)),
    if (gapend(1) < gapstart(1)),
        gapend = gapend(2:end);
    end;
    if (gapstart(end) > gapend(end)),
        gapstart = gapstart(1:end-1);
    end;
    
    for i = 1:length(gapstart),
        k = gapstart(i):gapend(i);
        cv(k) = cv(gapstart(i));
        
        ds(k) = ds(gapstart(i))/length(k);
    end;
end;
ds = ds(1:end-1);
       