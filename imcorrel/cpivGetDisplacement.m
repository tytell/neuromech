function [idx,idy, u,v, peak,err] = cpivGetDisplacement(corr, ind)

% only operate on indices ind.  If they don't give us ind, do all
% indices

px = size(corr,2);
py = size(corr,1);
ctrx = ceil(px/2 + 0.5);
ctry = ceil(py/2 + 0.5);

nx = size(corr,4);
ny = size(corr,3);
nphi = px*py;
npt = nx*ny;

if (nargin == 1),
	ind = 1:npt;
	isAll = 1;
else
	isAll = 0;
end;

corr = reshape(corr(:,:,ind), [nphi length(ind)]);
[peak, indcorr] = max(corr);

[idy,idx] = ind2sub([py px], indcorr);
bad = find(isnan(peak));
idx(bad) = NaN;
idy(bad) = NaN;

dxfrac = zeros(size(idx));
dyfrac = zeros(size(idy));
err = zeros(size(idx));
err(bad) = 16;

if (any(err == 0)),
	spx = [idx-1; idx; idx+1];
	k = find(all(spx >= 1) & all(spx <= px));
	spval = corr(sub2ind(size(corr), sub2ind([py px], repmat(idy(k), [3 1]), spx(:,k)), ...
					repmat(k, [3 1])));
					
	dxfrac(k) = cpivSubpixel(spval);
	lo = find(any(spx < 1));
	hi = find(any(spx > px));
	err(lo) = err(lo) + 1;
	err(hi) = err(hi) + 2;
	
	spy = [idy-1; idy; idy+1];
	k = find(all(spy >= 1) & all(spy <= py));
	spval = corr(sub2ind(size(corr), sub2ind([py px], spy(:,k), repmat(idx(k), [3 1])), ...
					repmat(k, [3 1])));
					
	dyfrac(k) = cpivSubpixel(spval);
	lo = find(any(spy < 1));
	hi = find(any(spy > py));
	err(lo) = err(lo) + 4;
	err(hi) = err(hi) + 8;
end;
	
u = idx + dxfrac - ctrx;
v = idy + dyfrac - ctry;
if (isAll)
	u = reshape(u, [ny nx]);
	v = reshape(v, [ny nx]);
	idx = reshape(idx, [ny nx]);
	idy = reshape(idy, [ny nx]);
	peak = reshape(peak, [ny nx]);
	err = reshape(err, [ny nx]);
end;

