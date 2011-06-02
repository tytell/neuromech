function [mxs,mys] = smoothEelMidline2(t,mx,my,eellen, serr,terr)
% function [mxs,mys] = smoothEelMidline2(t,mx,my,eellen, serr,terr)
% Smooths a midline using spaps over time and space simultaneously.  Checks
% the total length of the midline.  serr and terr are the errors along
% space and time, respectively.
%
% Mercurial revision hash: $Revision: 0570ba7d7b2b $ $Date: 2010/08/16 19:33:05 $
% Copyright (c) 2010, Eric Tytell

nfr = size(mx,2);
npt = size(mx,1);

s = [zeros(1,nfr); cumsum(sqrt(diff(mx).^2 + diff(my).^2))];

fprintf('Digitized length in pixels is %.2f+-%.2f = %.1f%% of input (%.2f)\n', ...
		nanmean(s(end,:)), nanstd(s(end,:))/sqrt(sum(isfinite(s(end,:)))), ...
		nanmean(s(end,:))/eellen*100, eellen);

if (nanstd(s(end,:))/eellen  > 0.05),
	warning('Lengths don''t match up (either digitized length varies or is different from input).');
end;

ds = diff(s);
ds0 = eellen/(npt-1);
if (any((ds - ds0)/ds0 > 0.05)),
	error('Distance along midline increases weirdly.');
end;

s = nanmedian(s,2);

k = find(all(isfinite(mx) & isfinite(my)));
XY = cat(1,shiftdim(mx(:,k),-1),shiftdim(my(:,k),-1));
sp = spaps({s,t(k)}, XY, {serr^2*range(s)*length(t(k)) terr^2*range(t(k))*length(s)});
XYs = fnval(sp,{s,t(k(1):k(end))});

mxs = NaN(size(mx));
mys = NaN(size(my));
mxs(:,k(1):k(end)) = squeeze(XYs(1,:,:));
mys(:,k(1):k(end)) = squeeze(XYs(2,:,:));

if (range(mxs(:)) > range(mys(:))),
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
plot(xx(:,q),yy(:,q),'.');
addplot(xxs(:,q),yys(:,q),'-');
xlabel('X position (mm)');
ylabel('Y position (mm)');
title('Three midlines at different times');

subplot(2,1,2);
q = [2 round(npt/2) npt-1];
plot(t,yy(q,:),'.');
addplot(t,yys(q,:),'k-');
xlabel('Time (s)');
ylabel('Y position (mm)');
title('Head, midbody and tail Y position over time');

