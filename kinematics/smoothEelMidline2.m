function [mxs,mys] = smoothEelMidline2(t,mx,my,eellen, serr,terr, varargin)
% function [mxs,mys] = smoothEelMidline2(t,mx,my,eellen, serr,terr)
% Smooths a midline using spaps over time and space simultaneously.  Checks
% the total length of the midline.  serr and terr are the errors along
% space and time, respectively.
%
% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2010, Eric Tytell

opt.discardbadspacing = true;
opt.offscreen = 'head';
opt = parsevarargin(opt, varargin, 7);

nfr = size(mx,2);
npt = size(mx,1);

s = [zeros(1,nfr); cumsum(sqrt(diff(mx).^2 + diff(my).^2))];

fprintf('Digitized length in pixels is %.2f+-%.2f = %.1f%% of input (%.2f)\n', ...
		nanmean2(s(end,:)), nanstd2(s(end,:))/sqrt(sum(isfinite(s(end,:)))), ...
		nanmean2(s(end,:))/eellen*100, eellen);

if (nanmedian2(s(end,:))/eellen  > 0.05),
	warning('Lengths don''t match up (either digitized length varies or is different from input).');
    actlen = nanmedian2(s(end,:));
else
    actlen = eellen;
end

ds = diff(s);
ds0 = actlen/(npt-1);
good = all(isfinite(mx) & isfinite(my));
goodspacing = all(abs((ds - ds0)/ds0) <= 0.2);

if (opt.discardbadspacing && any(~goodspacing))
    nbad = sum(good & ~goodspacing);
    warning('%d (%d%%) sequences discarded because of tracking errors', ...
        nbad, round(nbad/sum(good)*100));
    good = good & goodspacing;
end

if all(goodspacing)
    s = nanmedian2(s,2);

    k = find(good);
    XY = cat(1,shiftdim(mx,-1),shiftdim(my,-1));
    sp = spaps({s,t(k)}, XY(:,:,k), {serr^2*range2(s)*length(t(k)) terr^2*range2(t(k))*length(s)});
    XYs = NaN(size(XY));
    XYs(:,:,k(1):k(end)) = fnval(sp,{s,t(k(1):k(end))});
    offscreenind = NaN;
else
    switch opt.offscreen
        case 'head'
            s = s - repmat(s(end,:),[npt 1]) + actlen;
            isheadoff = true;
        case 'tail'
            isheadoff = false;
        otherwise
            error('Unrecognized offscreen option');
    end
    
    s0 = (0:ds0:actlen)';
    XY = cat(1,shiftdim(mx,-1),shiftdim(my,-1));
    XYs1 = NaN(size(XY));
    XYs = NaN(size(XY));
   
    len = range(s);
    
    offscreenind = first(abs(len - actlen) >= 0.5*ds0);
    
    %first spatial smoothing
    k = find(good);
    for ii = 1:length(k)
        i = k(ii);
        
        sp = spaps(s(:,i), XY(:,:,i), serr^2*range2(s(:,i)));
        
        if (abs(len(i) - actlen) < 0.5*ds0)
            s1 = s(:,i);
            srng = true(size(s1));
        else
            srng = (s0 >= min(s(:,i))) & (s0 <= max(s(:,i)));
            s1 = s0(srng);
            if ((s1(end) >= s(end,i)) && (s1(end) - s(end,i) < 0.5*ds0))
                s1(end) = s(end,i);
            elseif (length(s1) < size(XY,2))
                s1(end+1) = s(end,i);
                srng(length(s1)) = true;
            end
        end
        XYs1(:,srng,i) = fnval(sp,s1);
    end
    
    %then temporal smoothing
    for i = 1:size(XY,2)
        goodt = all(isfinite(XYs1(:,i,:)),1);
        k = find(goodt);
        
        sp = spaps(t(goodt),XYs1(:,i,goodt), terr^2*range2(t(goodt)));
        XYs(:,i,k(1):k(end)) = fnval(sp,t(k(1):k(end)));
    end
end
mxs = NaN(size(mx));
mys = NaN(size(my));
mxs = squeeze(XYs(1,:,:));
mys = squeeze(XYs(2,:,:));

if (range2(mxs(:)) > range2(mys(:))),
	ishoriz = 1;
	xx = mx;
	xxs = mxs;
	yy = my;
	yys = mys;
else
	ishoriz = 0;
	xx = -my;
	xxs = -mys;
	yy = mx;
	yys = mxs;
end;

subplot(2,1,1);
q = round([1 2 3]*nfr/4);
plot(xx(:,q),yy(:,q),'o');
addplot(xxs(:,q),yys(:,q),'-');
xlabel('X position (mm)');
ylabel('Y position (mm)');
title('Three midlines at different times');

subplot(2,1,2);
q = [2 round(npt/2) npt-1];
plot(t,yy(q,:),'o');
addplot(t,yys(q,:),'k-');

if (isfinite(offscreenind))
    vertplot(t(offscreenind),'k--');
    yl = ylim;
    text(t(offscreenind),yl(2), 'Off screen');
end

xlabel('Time (s)');
ylabel('Y position (mm)');
title('Head, midbody and tail Y position over time');

