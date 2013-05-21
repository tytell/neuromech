function [mx,my,good] = getMidlineExt(avi,frames, varargin)
% function [mx,my,good] = getMidlineExt(avi,frames, hx,hy,n, ...
%                                       width,len, maxang, ...
%                                       options)
% or
%
% function [mx,my,good] = getMidlineExt(avi,frames, hx,hy, tx,ty, n, ...
%                                       width,len, maxang, ...
%                                       options)
%
% If only hx and hy are passed, it attempts to progress up the body from
% head to tail, producing something of length len.
% If hx,hy and tx,ty are passed, it tries to subdivide the region between
% the two sets of points.
%
% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2010, Eric Tytell

hx = varargin{1};
hy = varargin{2};
if (all(size(varargin{3}) == size(hx)) && ...
        all(size(varargin{4}) == size(hx))),
    tx = varargin{3};
    ty = varargin{4};
    i = 5;
    method = 2;
else
    i = 3;
    tx = [];
    ty = [];
    method = 1;
end;
n = varargin{i};
width = varargin{i+1};
len = varargin{i+2};
maxang = varargin{i+3};

opts = varargin(i+4:end);

mincorrel = 0.7;
peakheight = 0.95;
maxpeakwidth = 1;

if (isnumeric(avi)),
    isavi = false;
    I = avi;
    if (isempty(frames)),
        frames = 1:size(I,3);
    end;
    nfr = size(I,3);
    vid = [];
else
    isavi = true;
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

invert = false;
eqhist = false;
quiet = false;
isvert = true;
debug = false;
enforcelen = true;
useEdges = false;
useHomogeneity = true;
headOutOfFrame = false(size(hx));
facedir = 1;
sumtransect = 0;
prevmx = [];
prevmy = [];
i = 1;
while (i <= length(opts)),
    switch lower(opts{i}),
        case 'invert',
            invert = true;
            if ((i < length(opts)) && islogical(opts{i+1})),
                invert = opts{i+1};
                i = i+1;
            end;
        case 'eqhist',
            eqhist = true;
            
        case 'quiet',
            quiet = true;
        case {'horizontal','horiz'},
            isvert = false;
        case {'facedown','faceright'},
            facedir = 1;
        case {'faceup','faceleft'},
            facedir = -1;
        case 'enforcetail',
            enforcelen = false;
        case 'debug',
            debug = 1;
        case 'sumtransect',
            sumtransect = opts{i+1};
            i = i+1;
        case 'mincorrel',
            mincorrel = opts{i+1};
            i = i+1;
        case 'useedges',
            useEdges = true;
        case 'noedges',
            useEdges = false;
        case 'usehomogeneity',
            useHomogeneity = true;
        case 'nohomogeneity',
            useHomogeneity = false;
        case 'headoutofframe',
            headOutOfFrame = opts{i+1};
            i = i+1;

        case 'previousmidline',
            prevmx = opts{i+1};
            prevmy = opts{i+2};
            i = i+2;
            
        otherwise,
            error('Unrecognized option %s',opts{i});
    end;
    i = i+1;
end;

if (any(size(headOutOfFrame) ~= size(hx))),
    error('Head out of frame data should be the same size as the head position');
end;
%if at any point the head leaves the frame, then we need method 3 - a
%combined one
if (any(headOutOfFrame)),
    method = 3;
end;

if (~quiet),
    set(gcf,'DoubleBuffer','on');
elseif (nfr > 1),
    timedWaitBar(0,'Finding midlines...');
end;

width = width.*len;
if (length(width) ~= n),
    width = spline(linspace(0,len,length(width)),width, ...
        linspace(0,len,n));
end;
width(width < 1) = 1;

seglen = len/(n-1);
perpd = 2*tan(maxang);
s = (0:seglen:len)';

actlen = [];
actseglen = [];

curvemax = maxang/(2*seglen);

%generate a first straight midline
if (~isempty(prevmx)),
    mx1 = prevmx;
    my1 = prevmy;
elseif (method == 2),
    mx1 = linspace(hx(frames(1)),tx(frames(1)),n)';
    my1 = linspace(hy(frames(1)),ty(frames(1)),n)';
elseif (isvert),
    mx1 = repmat(hx(frames(1)),[n 1]);
    my1 = hy(frames(1)) - facedir*s;
else
    mx1 = hx(frames(1)) - facedir*s;
    my1 = repmat(hy(frames(1)),[n 1]);
end;

mx = NaN([n nfr]);
my = NaN([n nfr]);
good = false(n,nfr);

%run through all the frames
for j = 1:length(frames),
    fr = frames(j);

    %skip frames with NaN head positions
    if (isnan(hx(fr)) || (~isempty(tx(fr)) && isnan(tx(fr)))),
        continue;
    end;

    %load the frame
    if (isavi),
        I1 = im2double(read(vid, fr));
        
        if (size(I1,3) > 1),
            I1 = I1(:,:,1);
        end;
    else
        I1 = I(:,:,fr);
    end;
    
    if (invert),
        I1 = 1-I1;
    end;
    if (eqhist),
        I1 = histeq(I1);
    end;
    
    %starting index and previous
    i = 2;
    p = 1;

    %move the previous frame's midline to the current head and/or tail
    %position
    switch method,
        case 1,
            %we don't have a tail position, so just move the head
            mx1 = mx1 - mx1(1) + hx(fr);
            my1 = my1 - my1(1) + hy(fr);

        case 2,
            %rotate and scale so that the previous midline fits between the
            %current head and tail positions
            [mx1,my1] = rotateAndScale(mx1,my1,1,hx(fr),hy(fr),n,tx(fr),ty(fr));

        case 3,
            if (~headOutOfFrame(fr)),
                [mx1,my1] = rotateAndScale(mx1,my1,1,hx(fr),hy(fr),n,tx(fr),ty(fr));
            else
                if (isempty(actlen)),
                    len = sum(sqrt(diff(mx(:,frames(1:j-1))).^2 + ...
                        diff(my(:,frames(1:j-1))).^2));
                    actlen = median(len);
                    actseglen = actlen / (n-1);
                end;
                
                onframe = isfinite(mx1);
                lenvisible = sum(sqrt(diff(mx1(onframe)).^2 + diff(my1(onframe)).^2));
                nvisible = floor(lenvisible / actseglen);
                i = n - nvisible + 1;
                p = i-1;

                mx1(1:i-2) = NaN;
                my1(1:i-2) = NaN;
                [mx1,my1] = rotateAndScale(mx1,my1, p,hx(fr),hy(fr), n,tx(fr),ty(fr));
%                 [mx1,my1] = rotateNotScale(mx1,my1,n,tx(fr),ty(fr),i,hx(fr),hy(fr));                
            end;
    end;

    p2 = NaN;
    good1 = false(n,1);
    good1(p) = 1;                        % head position is good
    prevseg = [NaN; NaN];

    %now step through the midline positions and try to find bright regions
    %that match the width
    done = 0;
    while ~done,
        dxds = (mx1(i)-mx1(p))/((i-p)*seglen);
        dyds = (my1(i)-my1(p))/((i-p)*seglen);

        %create a transect perpendicular to the midline, and proportional
        %to the number of steps we've gone since the last time we found a
        %good match
        step = perpd*(i-p)*seglen;
        perps = (-step-1.5*width(i)):(step+1.5*width(i));
        perps = perps + (step+1.5*width(i)-perps(end))/2;
        perpx = mx1(i) - perps*dyds;
        perpy = my1(i) + perps*dxds;

        trans = interp2(I1,perpx,perpy, '*nearest');

        %possibly create several transects nearby and sum them to try to
        %average out noise
        if (sumtransect > 0),
            for d = 1:sumtransect,
                perpx2 = mx1(i) + d*dxds - perps*dyds;
                perpy2 = my1(i) + d*dyds + perps*dxds;
                trans2 = interp2(I1,perpx2,perpy2,'*nearest');
                trans = trans+trans2;
                perpx2 = mx1(i) - d*dxds - perps*dyds;
                perpy2 = my1(i) - d*dyds + perps*dxds;
                trans2 = interp2(I1,perpx2,perpy2,'*nearest');
                trans = trans+trans2;
            end;
            trans = trans/(2*sumtransect + 1);
        end;

        %now generate a profile of the ideal image -- very bright and white
        %in the center, with the defined width, and black everywhere else
        profs = -0.5*width(i):0.5*width(i);
        profs = profs + (0.5*width(i)-profs(end))/2;
        profs = [(-3:-1)+profs(1) profs profs(end)+(1:3)];
        profile = zeros(size(profs));
        profile(abs(profs) < (width(i)+2)/2) = 0.5;
        profile(abs(profs) < width(i)/2) = 1;

        %convolve the ideal profile to determine the point of maximum
        %correlation
        C = conv(profile,trans)/width(i);

        off = 1:(length(perps)+length(profs)-1);
        off = off - (off(end)-1)/2;

        %grab the center of the correlation function
        doff = off(2)-off(1);
        a = (perps(1)-off(1))/doff;
        if (a - floor(a) == 0),
            k = a+(0:length(perps)-1);
            C = C(k);
            off = off(k);
        elseif (a - floor(a) == 0.5),
            k = floor(a)+(0:length(perps)-1);
            C = (C(k)+C(k+1))/2;
            off = (off(k)+off(k+1))/2;
        end;

        w2 = ceil(0.9*width(i)/2);
        if (w2 == 0),
            %can't have zero width regions
            w2 = 1;
        end;

        %find the maximum correlation
        [maxC,ind] = max(C);

        %if reflections have maxed out the image, we may end up with a
        %rather flat peak, so look for where the peak drops off on either
        %side and take the center between the two drop offs, not
        %(necessarily) the absolute maximum
        a = ind;
        amin = w2;
        while ((a > amin) && (C(a) > peakheight*maxC)),
            a = a-1;
        end;
        b = ind;
        bmax = length(perps)-w2;
        while ((b < bmax) && (C(b) > peakheight*maxC)),
            b = b+1;
        end;
        ctr = round((a+b)/2);
        if (useHomogeneity && (w2 > 1)),
            %also check for relatively homogeneous areas with slightly less than
            %the appropriate width
            transgrad = deriv(perps,trans);
            rms = repmat(NaN,size(trans));

            for rind = w2:length(transgrad)-w2,
                k = (-w2+1:w2) + rind;
                rms(rind) = sqrt(sum(transgrad(k).^2)/width(i));
            end;
            if (any(rms ~= 0)),
                rms = rms/nanmedian(rms);
                if (any(rms == 0)),
                    rms((rms == 0) & (trans == 1)) = min(rms(rms > 0));
                    rms((rms == 0) & (trans < 1)) = Inf;
                end;

                %maximize correlation and homogeneity (= 1/rms)
                [corrhomog,newctr] = max(C ./ rms);
                if ((C(newctr) > 0.2) && (C(newctr) > 0.75*maxC)),
                    maxC = corrhomog;
                    ctr = newctr;
                end;
            end;
        end;

        if (useEdges),
            %also check for edges (represented by the gradient)
            transgrad = deriv(perps,trans);

            %find strong changes in intensity
            up = makecol(find(transgrad > 0.1));
            down = makerow(find(transgrad < -0.1));

            %that are separated by the defined width +- 10%
            dist = repmat(down,[length(up) 1]) - ...
                repmat(up,[1 length(down)]);
            [l,r] = find(abs((dist-width(i))/width(i)) < 0.1);
            l = up(l);
            r = down(r);

            %also grab the edges that are closest to our current segment
            [q,l0] = min(abs(ctr - width(i)/2 - up));
            [q,r0] = min(abs(ctr + width(i)/2 - down));
            l0 = up(l0);
            r0 = down(r0);

            %if we find good edges, then check whether there's a correlation
            %peak between them
            if (~isempty(l)),
                weight = transgrad(l)-transgrad(r);

                [w,edgectr] = max(weight);

                %check whether the height of the correlation peak plus the
                %strength of the edges is greater than the original correlation
                %peak plus its edges
                if (w + C(edgectr) >= transgrad(l0) - transgrad(r0) + maxC),
                    ctr = edgectr;
                    maxC = w + C(ctr);

                    a = ctr;
                    b = ctr;
                end;
            end;
        end;

        %calculate the angular change
        curseg = [(mx1(i)-off(ctr)*dyds-mx1(i-p)) ...
            (my1(i)+off(ctr)*dxds-my1(i-p))];
        curseg = curseg./sqrt(sum(curseg.^2));
        angdiff = acos(curseg*prevseg)/(i-p2-1);

        if ((maxC >= mincorrel) && (a > amin) && (b < bmax) && ...
                (isnan(angdiff) || (j == 1) || (angdiff < maxang))),
            opp = off(ctr);
            adj = (i-p)*seglen;

            hyp = sqrt(opp^2 + adj^2);

            if (i < n),
                switch method,
                    case 1,
                        rotx = adj/hyp*(mx1(i:end)-mx1(p)) - ...
                            opp/hyp*(my1(i:end)-my1(p)) ...
                            + mx1(p);
                        roty = opp/hyp*(mx1(i:end)-mx1(p)) + ...
                            adj/hyp*(my1(i:end)-my1(p)) ...
                            + my1(p);

                        mx1(i:end) = rotx;
                        my1(i:end) = roty;

                    case {2,3},
                        rotx = adj/hyp*(mx1(i)-mx1(p)) - ...
                            opp/hyp*(my1(i)-my1(p)) ...
                            + mx1(p);
                        roty = opp/hyp*(mx1(i)-mx1(p)) + ...
                            adj/hyp*(my1(i)-my1(p)) ...
                            + my1(p);

                        [mx1,my1] = ...
                            rotateAndScale(mx1,my1, i,rotx,roty, n,tx(fr),ty(fr));
                end;
            else
                rotx = adj/hyp*(mx1(i)-mx1(p)) - ...
                    opp/hyp*(my1(i)-my1(p)) ...
                    + mx1(p);
                roty = opp/hyp*(mx1(i)-mx1(p)) + ...
                    adj/hyp*(my1(i)-my1(p)) ...
                    + my1(p);

                mx1(i) = rotx;
                my1(i) = roty;
            end;

            p2 = p;
            p = i;
            prevseg = curseg';

            good1(i) = 1;
        end;

        if (i == n),
            done = 1;
        else
            i = i+1;
        end;
    end;

    if (~enforcelen || (~isempty(tx) && ~good1(end))),
        mx1(end) = tx(fr);
        my1(end) = ty(fr);
        good1(end) = 1;
    end;

    s1 = [0; cumsum(sqrt(diff(mx1(good1)).^2 + diff(my1(good1)).^2))];

    if (headOutOfFrame(fr)),
        %we only want to extrapolate for those points that are actually on
        %the frame, so figure out how many of those there are
        seglen1 = sqrt(diff(mx1(good1)).^2 + diff(my1(good1)).^2);
        sfromtail = [0; cumsum(seglen1(end:-1:1))];
        seq = 0:actseglen:sfromtail(end);
        
        sfromtail = actlen - sfromtail(end:-1:1);
        s1 = sfromtail;
        s2 = repmat(NaN,[n 1]);
        s2(n-length(seq)+1:end) = actlen - seq(end:-1:1);
    elseif (enforcelen)
        s2 = s;
    else
        s2 = linspace(0,s1(end),n);
    end;

    sp = csape(s1,[mx1(good1) my1(good1)]', 'variational');
    mxy = fnval(sp,s2(isfinite(s2)));
    mx1 = repmat(NaN,[n 1]);
    my1 = repmat(NaN,[n 1]);
    mx1(isfinite(s2)) = mxy(1,:)';
    my1(isfinite(s2)) = mxy(2,:)';

    ds = sqrt(diff(mx1).^2 + diff(my1).^2);

    if (~good1(end)),
        if (all(good1(end-3:end-1))),
            dxyds = fnval(fnder(sp),s(end-1));
            dxyds2 = fnval(fnder(sp,2),s(end-1));

            mx1(end) = mx1(end-1) + seglen*dxyds(1) + ...
                0.5*seglen^2*dxyds2(1);
            my1(end) = my1(end-1) + seglen*dxyds(2) + ...
                0.5*seglen^2*dxyds2(2);
        elseif (all(good1(end-4:end-2))),
            dxyds = fnval(fnder(sp),s(end-2));
            dxyds2 = fnval(fnder(sp,2),s(end-2));

            mx1(end-1:end) = mx1(end-2) + (1:2)*seglen*dxyds(1) + ...
                0.5*((1:2)*seglen).^2*dxyds2(1);
            my1(end-1:end) = my1(end-2) + (1:2)*seglen*dxyds(2) + ...
                0.5*((1:2)*seglen).^2*dxyds2(2);
        end;
    end;

    if (~quiet),
        imshow6(I1,'n');
        drawx = mx1;
        drawy = my1;
        drawx(~good1) = NaN;
        drawy(~good1) = NaN;
        addplot(drawx,drawy,'r.-',mx1(~good1),my1(~good1),'g.');
        if (method == 2),
            addplot(tx(fr),ty(fr),'ro');
        end;
        drawnow;

        set(gcf,'Name',sprintf('Frame %d/%d (%d%%)',fr,nfr, ...
            round(j/length(frames)*100)));
    elseif (nfr > 1),
        timedWaitBar(j/length(frames));
    end;

    mx(:,fr) = mx1;
    my(:,fr) = my1;
    good(:,fr) = good1;
end;

if (quiet && (nfr > 1)),
    timedWaitBar(1);
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
function [x,y] = rotateNotScale(x,y, a,xa,ya, edge,xedge,yedge)

%vector of old midline
dx1 = x(edge)-x(a);
dy1 = y(edge)-y(a);
%vector to rotate to
dx2 = xedge-xa;
dy2 = yedge-ya;
%distances for normalization
olddist = sqrt(dx1^2 + dy1^2);
newdist = sqrt(dx2^2 + dy2^2);

%take the dot product
dotprod = (dx1*dx2 + dy1*dy2) / (olddist*newdist);
%keeping track of which direction we need to rotate
sgn = sign(dx1*dy2 - dx2*dy1);

if (dotprod > 1),                       % rounding error problem
    dotprod = 1;
end;

ang = sgn*acos(dotprod);

x2 = ((x - x(a))*cos(ang) - (y - y(a))*sin(ang)) + xa;
y2 = ((x - x(a))*sin(ang) + (y - y(a))*cos(ang)) + ya;

x = x2;
y = y2;

