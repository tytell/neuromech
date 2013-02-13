function [mx,my] = getImageMidline(I, hx,hy, tx,ty, n, width)

maxang = pi/4;

w2 = round(width/2);
profile = [zeros(w2,1); ones(width,1); zeros(w2,1)];

ndiv = ceil(log2(n));
mindiv = 1;

% initial boundaries are just the head and tail points
mx = [hx; tx];
my = [hy; ty];

% distance between points
a = sqrt(diff(mx).^2 + diff(my).^2);
len = a;

for j = 1:ndiv,
	nx = length(mx) - 1;
	len = len/2;
	
	% angle between each segment
	theta = atan2(diff(my),diff(mx));
	a = a / 2;

	% points along the perpendicular bisector of each segment
	% perplen is the approximate number of pixels in the bisector
	% new points are x = x0 + l cos(alpha - theta) and y = y0 + l sin(alpha - theta)
	% old points (x0,y0) are in mx,my as a column vector
	% alpha is a row vector
	perplen = round(a*tan(maxang));
	alpha = linspace(-maxang/2, maxang/2, perplen);
	newx = repmat(mx(1:end-1),[1 perplen]) + a./repmat(cos(alpha), [nx 1]) .* ...
			cos(repmat(alpha, [nx 1]) + repmat(theta, [1 perplen]));
	newy = repmat(my(1:end-1),[1 perplen]) + a./repmat(cos(alpha), [nx 1]) .* ...
			sin(repmat(alpha, [nx 1]) + repmat(theta, [1 perplen]));
	ind = repmat(1:perplen, [nx 1]);
	
	% interpolate intensities along the bisectors
	% and choose the maximum as the new point
	newI = interp2(I, newx, newy, '*linear');
	
	off = (perplen - length(profile))/2;
	cx = -ceil(off):floor(off);

	for k = 1:nx,
		C = normxcorr1(profile(k,:), newI(k,:));
		[maxC, ind] = max(C);
		
		newx(k,1) = newx(k,ceil(perplen/2) + cx(ind));
		newy(k,1) = newy(k,ceil(perplen/2) + cx(ind));
	end;
	
	newx = newx(:,1)';
	newy = newy(:,1)';
	
	% merge new points with the old ones.  NB: to be ordered correctly
	% the merged vector should be [mx(1) newx(1) my(2) newy(2), etc.]
	% we do this my flattening a two row matrix
	newx = [mx(1:end-1) newx]';
	newy = [my(1:end-1) newy]';
	mx = [newx(:); mx(end)];
	my = [newy(:); my(end)];
	nx = length(mx) - 1;
	
	profile = newI;
end;

