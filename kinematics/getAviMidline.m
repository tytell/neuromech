function [mx,my] = getImageMidline(I, hx,hy, tx,ty, n)

w2 = round(width/2);
profile = [zeros(1,w2) ones(1,width) zeros(1,w2)];

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
	perplen = round(a*tan(maxang))*2;
	alpha = linspace(-maxang/2, maxang/2, perplen);
	newx = repmat(mx(1:end-1),[1 perplen]) + a./repmat(cos(alpha), [nx 1]) .* ...
			cos(repmat(alpha, [nx 1]) + repmat(theta, [1 perplen]));
	newy = repmat(my(1:end-1),[1 perplen]) + a./repmat(cos(alpha), [nx 1]) .* ...
			sin(repmat(alpha, [nx 1]) + repmat(theta, [1 perplen]));
	ind = repmat(1:perplen, [nx 1]);
	
	% interpolate intensities along the bisectors
	% and choose the maximum as the new point
	newI = interp2(I, newx, newy, '*linear');
	
	
	[maxI,pos] = max(newI, [], 2);
	ind(newI < repmat(maxI,[1 perplen])) = NaN;
	pos = round(nanmean(ind,2));
	
	newx = newx(sub2ind(size(newx), (1:size(newx,1))', pos));
	newy = newy(sub2ind(size(newy), (1:size(newy,1))', pos));
	
	% merge new points with the old ones.  NB: to be ordered correctly
	% the merged vector should be [mx(1) newx(1) my(2) newy(2), etc.]
	% we do this my flattening a two row matrix
	newx = [mx(1:end-1) newx]';
	newy = [my(1:end-1) newy]';
	mx = [newx(:); mx(end)];
	my = [newy(:); my(end)];
	
	% now calculate the coordinates of each pixel along each segment
	sint = (0:round(len)-1)';
	ns = length(sint);
	nx = length(mx) - 1;
	mxall = repmat(mx(1:end-1)', [ns 1]) + repmat(sint, [1 nx]) .* repmat(diff(mx)'./a, [ns 1]);
	myall = repmat(my(1:end-1)', [ns 1]) + repmat(sint, [1 nx]) .* repmat(diff(my)'./a, [ns 1]);
	mxall = mxall(:);
	myall = myall(:);
	
	% and integrate the intensities along the line
	s = [0; cumsum(sqrt(diff(mxall).^2 + diff(myall).^2))];
	newint = trapz(s, interp2(I, mxall,myall, '*linear'));
end;

