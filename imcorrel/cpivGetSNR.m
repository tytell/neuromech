function [snr,rmsnoise,peaknoise] = cpivGetSNR(phi, idx,idy, sig, ind)

nx = size(phi,4);
ny = size(phi,3);

if (nargin == 5),
	npt = length(ind);
	isAll = 0;
else
	npt = nx*ny;
	ind = 1:npt;
	isAll = 1;
end;

px = size(phi,2);
py = size(phi,1);
phi1 = reshape(phi(:,:,ind), [px*py npt]);

% make them row vectors
ind = ind(:)';
if (prod(size(idx)) ~= npt),
	idx = idx(ind)';
	idy = idy(ind)';
else
	idx = idx(:)';
	idy = idy(:)';
end;

[peakx,peaky] = meshgrid(-1:1);
npk = 9;

idx = repmat(idx,[npk 1]) + repmat(peakx(:),[1 npt]);
idy = repmat(idy,[npk 1]) + repmat(peaky(:),[1 npt]);
ipt = repmat(1:npt,[npk 1]);

k = find(all(idx >= 1) & all(idx <= px) & all(idy >= 1) & all(idy <= py) & ...
			any(isfinite(phi1)));

if (~isempty(k)),
	peakind = sub2ind([py px npt], idy(:,k),idx(:,k),ipt(:,k));
	
	% add up just the peak
	% sig = repmat(NaN, [npt 1]);
	% sig(k) = sqrt(sum(phi2(peakind)) / 9);
	
	% now zero out the peak
	phi1(peakind) = 0;
	
	% and add up the rest
	rmsnoise = repmat(NaN, size(sig));
	rmsnoise(k) = nansum(phi1(:,k).^2) / (sum(isfinite(phi1(:,k))) - 9);
	rmsnoise = sqrt(rmsnoise);
	
	% find the peak noise value
	peaknoise = repmat(NaN, size(sig));
	peaknoise(k) = max(phi1(:,k));
	
	q = find(peaknoise ~= 0);
	snr = repmat(NaN, size(sig));
	snr(q) = sig(q)./peaknoise(q);
end;
