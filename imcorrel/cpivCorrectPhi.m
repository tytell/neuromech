function [phi, iscorrind, harterr] = cpivCorrectPhi(phi, u,v, snr, snrmin)

trydx = [0 0  1 -1]';
trydy = [1 -1 0  0]';

nx = size(u,2);
ny = size(u,1);
[x,y] = meshgrid(1:nx, 1:ny);

harterr = zeros(size(snr));

needscorr = find(isfinite(snr) & (snr < snrmin));
harterr(needscorr) = 1;

if (~isempty(needscorr)),
	x = x(needscorr)';
	y = y(needscorr)';
	npt = length(needscorr);
	
	corrx = repmat(x, [4 1]) + repmat(trydx, [1 npt]);
	corry = repmat(y, [4 1]) + repmat(trydy, [1 npt]);
	
	nearsnr = repmat(NaN, [4 npt]);
	
	k = find((corrx >= 1) & (corrx <= nx) & (corry >= 1) & (corry <= ny));
	nearsnr(k) = snr(sub2ind([ny nx], corry(k), corrx(k)));
	
	[maxnearsnr, corrdir] = nanmax(nearsnr);
	k = find(isfinite(corrdir) & isfinite(maxnearsnr));
	iscorrind = needscorr(k);
	harterr(iscorrind) = 0;
	
	corrx = x(k) + trydx(corrdir(k))';
	corry = y(k) + trydy(corrdir(k))';
	
	phi(:,:,iscorrind) = phi(:,:,iscorrind) .* phi(:,:, sub2ind([ny nx],corry,corrx));
else
	iscorrind = [];
end;

