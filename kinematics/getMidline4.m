function [mx,my,good] = getMidline4(avi,frames, hx,hy,tx,ty, mx1,my1, n,len,maxang, varargin)

opt.invert = false;
opt.quiet = false;
opt.subtractbackground = [];
opt.eqhist = false;
opt.ntransect = 5;
opt.gridsize = 10;

opt = parsevarargin(opt,varargin, 10);

if (isnumeric(avi)),
    error('Cannot take image matrices');
else
    vid = VideoReader2(avi);
    if (isempty(frames)),
        frames = 1:vid.NumberOfFrames;
    elseif (numel(frames) == 1),
        frames = 1:frames:vid.NumberOfFrames;
    end;
    nfr = vid.NumberOfFrames;
end;

if ((length(hx) < frames(end)) || any(~isfinite(hx(frames)))),
    error('Need head position for each frame');
end;

if (frames(1) == 1)
    frames = frames(2:end);
end
if (frames(end) == nfr)
    frames = frames(1:end-1);
end
iseven = all(diff(frames) == 1);

seglen = len/(n-1);
perpd = ceil(0.1*len);
maxoff = ceil(tan(maxang)*seglen);
s = (0:seglen:len)';

actlen = [];
actseglen = [];

curvemax = maxang/(2*seglen);

mx = NaN([n nfr]);
my = NaN([n nfr]);
good = false(n,nfr);

fr = frames(1)-1;
I1 = im2double(read(vid, fr));
I1 = processFrame(I1, opt);

fig = gcf;
clf;
hax = gca;
him = imshow(I1);

if (isempty(mx1))
    zoom on;
    input('Zoom to fish and press return');
    zoom off;
    addplot(hx(fr),hy(fr),'g+', tx(fr),ty(fr),'g+');
    fprintf('Click first midline (excluding head and tail)\n');
    [mx0,my0] = ginputb;
    
    mx0 = [hx(fr); mx0; tx(fr)];
    my0 = [hy(fr); my0; ty(fr)];
    
    s0 = [0; cumsum(sqrt(diff(mx0).^2 + diff(my0).^2))];
    
    sp = spaps(s0,[mx0 my0]', 0.5^2);
    
    mxy1 = fnval(sp, s);
    mx1 = mxy1(1,:)';
    my1 = mxy1(2,:)';
end;

hln = addplot(mx1,my1,'g+-');
hln(2) = addplot(mx1,my1,'bo-');
lim = [min(mx1)-0.5*len max(mx1)+0.5*len min(my1)-0.5*len max(my1)+0.5*len];
axis(hax,lim);

ntrans1 = floor(opt.ntransect/2);

for i = 1:length(frames)
    fr = frames(i);
    
    %load the frame
    if (iseven)
        Iprev = I1;
        I1 = im2double(read(vid, fr));
        I1 = processFrame(I1,opt);
    else
        I1 = im2double(read(vid, fr));
        I1 = processFrame(I1,opt);
    end
    
    %get the next frame
    %Inext = im2double(read(vid, fr+1));
    %Inext = processFrame(Inext,opt);
    
    %rotate and scale so that the previous midline fits between the
    %current head and tail positions
    [mx2,my2] = rotateAndScale(mx1,my1,1,hx(fr),hy(fr),n,tx(fr),ty(fr));

    dxds = zeros(n,1);
    dyds = zeros(n,1);
    dxds(2:end-1) = (mx1(3:end) - mx1(1:end-2))/(2*seglen);
    dyds(2:end-1) = (my1(3:end) - my1(1:end-2))/(2*seglen);
    mag = sqrt(dxds.^2 + dyds.^2);
    dxds(2:end-1) = dxds(2:end-1)./mag(2:end-1);
    dyds(2:end-1) = dyds(2:end-1)./mag(2:end-1);
    
    offx1 = mx2 - mx1;
    offy1 = my2 - my1;
    
    offpar = offx1.*dxds + offy1.*dyds;
    offperp = zeros(n,1);
    
    maxcorr = zeros(n,1);
    for pt = 2:n-1
        %create a transect perpendicular to the midline
        step = perpd; %*seglen;
        perps = -step:step;
        perps = perps - (perps(1)+perps(end))/2;        % center on 0
        perps = repmat(perps,[2*ntrans1+1 1]);
        perpn = (-ntrans1:ntrans1)';
        perpn = repmat(perpn,[1 length(perps)]);
        
        perpx = mx1(pt) - perps*dyds(pt) + perpn*dxds(pt);
        perpy = my1(pt) + perps*dxds(pt) + perpn*dyds(pt);

        trans1 = interp2(I1,perpx,perpy, '*linear');
        transprev = interp2(Iprev,perpx,perpy, '*linear');
        %transnext = interp2(Inext,perpx,perpy, '*linear');
        
        C1 = normxcorr2(transprev,trans1);
        k = (-maxoff:maxoff)+size(trans1,2);
        [maxC1, imax] = max(flatten(C1(:,k)));
        [npeak1,speak1] = ind2sub([size(C1,1) length(k)], imax);
        speak1 = speak1+k(1)-1;
        frac1 = subpixel(C1(npeak1,speak1+(-1:1))');
        speak1 = speak1 - size(trans1,2) + frac1;
        
%         C2 = normxcorr2(trans1,transnext);
%         [maxC2, imax] = max(flatten(C2(:,k)));
%         [npeak2,speak2] = ind2sub([size(C2,1) length(k)], imax);
%         speak2 = speak2+k(1)-1;
%         frac2 = subpixel(C2(npeak2,speak2+(-1:1))');
%         speak2 = speak2 - size(trans1,2) + frac2;

        offperp(pt) = speak1;
        maxcorr(pt) = maxC1;
    end
    
    mx2 = mx1 + offpar.*dxds - offperp.*dyds;
    my2 = my1 + offpar.*dyds + offperp.*dxds;
    mx2([1 end]) = [hx(fr); tx(fr)];
    my2([1 end]) = [hy(fr); ty(fr)];
    
    ang1 = atan2(diff(my2), diff(mx2));
    curve1 = [0; diff(ang1) / (2*seglen); 0];
    dcurveds1 = (abs(curve1(3:end) - curve1(2:end-1)) + ...
        abs(curve1(2:end-1) - curve1(1:end-2)))/(2*seglen);
    
    weight1 = 1./dcurveds1;
    weight1(weight1 > 1e6) = 1e6;
    weight1 = [2*max(weight1); weight1; 2*max(weight1)];
    weight1 = weight1 / sum(weight1);
    
    s1 = [0; cumsum(sqrt(diff(mx2).^2 + diff(my2).^2))];
    sp = spaps(s1', [mx2 my2]', -0.8, weight1);
    
    mxy = fnval(sp, s);
    mx1 = mxy(1,:)';
    my1 = mxy(2,:)';

    mx1(end) = tx(fr);
    my1(end) = ty(fr);
    
    mx(:,fr) = mx1;
    my(:,fr) = my1;
    
    if (~opt.quiet)
        set(him,'CData',I1);
        set(hln(1),'XData',mx(:,fr-1),'YData',my(:,fr-1));
        set(hln(2),'XData',mx1,'YData',my1);

        lim = [min(mx1)-0.5*len max(mx1)+0.5*len min(my1)-0.5*len max(my1)+0.5*len];
        axis(hax,lim);
        drawnow;
    end

end


function I = processFrame(I, opt)

if (size(I,3) > 1),
    I = I(:,:,1);
end;
if (~isempty(opt.subtractbackground))
    I = I - opt.subtractbackground;
    I = I - min(I(:));
end

if (opt.invert),
    I = 1-I;
end;
if (opt.eqhist)
    I = histeq(I);
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

% -------------------------------------------------------
function frac = subpixel(val)

% remove negative values
m = min(val);
isneg = m <= 0;
if (any(isneg))
	val(:,isneg) = val(:,isneg) - repmat(m(isneg), [3 1]) + 0.0001;
end;

lval = log(val);

frac = (lval(1,:) - lval(3,:))./(2*(lval(1,:) + lval(3,:) - 2*lval(2,:)));

% -------------------------------------------------------
function showCorrelation(transprev,trans1, C1, speak1,maxoff)

x = (1:size(trans1,2)) - size(trans1,2)/2;
y = (1:size(trans1,1)) - size(trans1,1)/2;

subplot(2,2,1);
imagesc(x,y,transprev);
addplot(0,0,'r+','MarkerSize',12);

subplot(2,2,3);
imagesc(x,y,trans1);
addplot([0 speak1],[0 0],'r+','MarkerSize',12);

lag = (1:size(C1,2)) - size(trans1,2);

subplot(4,2,[4 6]);
plot(lag,C1');
vertplot(speak1,'k-');
xlim([-maxoff maxoff]);



    
    