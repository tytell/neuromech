function [mx,my, lx,ly, rx,ry] = getMidline(I, hx,hy, tx,ty, n,width, maxang)
% function [mx,my, lx,ly, rx,ry] = getMidline(I, hx,hy, tx,ty, n,width, [maxang])
%
% Find fish midline between a head and a tail point.  Progressively
% subdivides the region between the head and the tail into halves and
% looks along the perpendicular bisector for bright regions.
% Uses a 1D cross correlation along the bisector to find the fish.  It
% assumes the fish has a square intensity profile, width pixels wide.
% The cross-correlation peak between the square intensity profile along
% the perpendicular bisector is assumed to be the midline position.
% Then it bisects each segment again and repeats the process, until it
% has found n points.
%
% Parameters:
%	I - Image matrix (double, not uint8)
%	hx,hy - Head position in pixels
%	tx,ty - Tail position in pixels
%	n - Number of points along the midline to find
%	width - Width of the fish in fish lengths (not pixels).  Must contain
%		evenly spaced points from the head to the tail, but it can have
%		any number of points.
%	maxang - Optional argument giving the maximum angle in radians between
%		any two segments.  It can reduce processing time if it's low.

centerBiasStrength = [];

% Find the number of subdivisions we need to make to get the right
% number of points
ndiv = ceil(log2(n));
npt = 2^ndiv + 1;

if ((nargin == 7) | isempty(maxang)),
	maxang = 75*pi/180;
end;

% The midline (mx,my) starts with just two points, the head and the tail
% Also calculate the total "length" between the points
mx = [hx; tx];
my = [hy; ty];
total = sqrt(diff(mx).^2 + diff(my).^2);

if (isempty(centerBiasStrength)),
    centerBiasStrength = 1;
end;

% reduce the image size so that processing is quicker

% First calculate the extent of the region we'll need to look at
% The first bisection results in these points:
angx = mean(mx) - [-1;1]*tan(maxang/2)*diff(my)/2;
angy = mean(my) + [-1;1]*tan(maxang/2)*diff(mx)/2;

% Find the maximum range
w = range([mx; angx]) * 1.1;
h = range([my; angy]) * 1.1;
ctrx = (max([mx; angx]) + min([mx; angx]))/2;
ctry = (max([my; angy]) + min([my; angy]))/2;

% And define the region that covers that range
rng = round([ctrx-w/2 ctrx+w/2 ctry-h/2 ctry+h/2]);
rng(rng < 1) = 1;
if (rng(2) > size(I,2)),
	rng(2) = size(I,2);
end;
if (rng(4) > size(I,1)),
	rng(4) = size(I,1);
end;

% reduce the image size
I = I(rng(3):rng(4), rng(1):rng(2));
% and make sure the midline is in the new coordinates
mx = mx - rng(1)+1;
my = my - rng(3)+1;

% Interpolate width so it has the number of points we're going to calculate,
% and multiply it by the length so that it's in pixels
width = shiftdim(width)';
nw = length(width);
width = csapi(linspace(0,1,nw), width, 0:1/2^ndiv:1) * total;
maxw = max(width);
maxw2 = floor((maxw-1)/2);
% Center of the intensity profile
% NB: having the profile centered in the variable profile is very
% important, because if it's off center, then the cross correlation
% peak will be offset, too
prof0 = 2*maxw2 + 5;
% Distance along the profile
profdist = -(2*maxw2 + 3):(2*maxw2 + 3);
nprof = length(profdist);

% Build the intensity profile for each width
profile = zeros(length(width), nprof);
for i = 1:length(width),
	% integer width
	w2 = floor((width(i)-1)/2);
	% fractional width corresponds to lower intensity
	wfrac = width(i)-(2*w2 + 1);
	
	profile(i,prof0 + (-w2:w2)) = 1;
	profile(i,prof0 + [-w2-1 w2+1]) = wfrac/2;
end;

minlen = 1.5*nprof;

% Now do the subdivisions
for i = 1:ndiv,
	% Length of each segment
	len = sqrt(diff(mx).^2 + diff(my).^2);
	% And the normalized dx and dy of each segment (dx^2 + dy^2 = 1)
	dx = diff(mx)./len;
	dy = diff(my)./len;
	
	% Length of the perpendicular bisector
	perplen = mean(len)*tan(maxang/2);
	% It can't be less than minlen
	perplen(perplen < minlen) = minlen;
	
	% Actual distance along the bisector
	perpdist = -round(perplen/2):round(perplen/2);
	% Number of points
	nperp = size(perpdist,2);

	% Midpoints of each segment
	newx = mx(2:end) - dx.*len/2;
	newy = my(2:end) - dy.*len/2;
	nnew = length(newx);
	
	% And the coordinates of the bisectors
	perpx = repmat(newx, [1 nperp]) - ...
				repmat(perpdist, [nnew 1]).*repmat(dy, [1 nperp]);
	perpy = repmat(newy, [1 nperp]) + ...
				repmat(perpdist, [nnew 1]).*repmat(dx, [1 nperp]);
	
	% Use interp2 to interpolate image intensities along the bisectors
	xsect = interp2(I,perpx,perpy);
	% off defines the offset from the center of the bisector -- it'll be
	% used with the cross correlation
	off = (length(perpdist) - size(profile,2))/2;
	off = -ceil(off):floor(off);

	% Find the appropriate profile for each point
	profind = 1:(npt-1)/2^i:npt;
	profind = profind(2:2:end);
	
	% And run the cross correlation
	for j = 1:length(newx),
		C = xcorr1(profile(profind(j),:), xsect(j,:));
                bias = exp(-off.^2/max(off).^2*centerBiasStrength.^2);
                C = C.*bias;

		[maxcorr, ind] = max(C);
		
		% The maximum correlation value corresponds to a distance off(ind)
		% from the center of the bisector.  Set the new value to that
		% position
		newx(j) = newx(j) - dy(j)*off(ind);
		newy(j) = newy(j) + dx(j)*off(ind);
	end;
	
	% Integrate the new points into the midline
	mx2(1:2:2*length(mx)-1,1) = mx;
	mx2(2:2:2*length(mx)-1,1) = newx;
	my2(1:2:2*length(my)-1,1) = my;
	my2(2:2:2*length(my)-1,1) = newy;
	
	mx = mx2;
	my = my2;
end;

% Now we have a power of 2 points.  Interpolate them down to get
% only n points.
s0 = [0; cumsum(sqrt(diff(mx).^2 + diff(my).^2))];
s = linspace(s0(1),s0(end),n)';
mxy = csapi(s0', [mx my]', s');
mx = mxy(1,:)';
my = mxy(2,:)';

% Repeat the cross correlation process to find the edges of the fish
% at each of the interpolated midline points.  This time we only use
% half of the profile to find edges

% Normalized x and y distances between segments
dx = (mx(3:end) - mx(1:end-2))./(s(3:end) - s(1:end-2));
dy = (my(3:end) - my(1:end-2))./(s(3:end) - s(1:end-2));

perplen = minlen;

% Distance along the bisector
perpdist = -round(perplen/2):round(perplen/2);
nperp = size(perpdist,2);

% Bisector points
perpx = repmat(mx(2:end-1), [1 nperp]) - ...
			repmat(perpdist, [n-2 1]).*repmat(dy, [1 nperp]);
perpy = repmat(my(2:end-1), [1 nperp]) + ...
			repmat(perpdist, [n-2 1]).*repmat(dx, [1 nperp]);

% Build the half profiles.  Again, be careful that the point we're
% interested in (the edge) is at the center of the profile, and that
% the profile has the same number of points on either side of the
% edge
prof2 = repmat(NaN, [npt 2*maxw2 + 5]);
prof0 = maxw2+3;
prof2dist = -maxw2-2:maxw2+2;

for i = 1:npt,
	w2 = floor((width(i)-1)/2);
	wfrac = (width(i)-(2*w2 + 1))/2;

	prof2(i, prof0 + (-w2-2:-2)) = 1;
	prof2(i, prof0 - 1) = wfrac;
	prof2(i, prof0 + (0:w2+2)) = 0;
end;

% Interpolate image intensities
xsect = interp2(I,perpx,perpy);

% Offsets from the center
off = (length(perpdist) - size(prof2,2))/2;
off = -ceil(off):floor(off);

% Left and right points at the head and the tail are just the
% same as the midline
lx([1 n],1) = mx([1 n]);
rx([1 n],1) = mx([1 n]);
ly([1 n],1) = my([1 n]);
ry([1 n],1) = my([1 n]);

% Find the profiles that correspond to the new (non-power-of-2) n points
profind = round((1:n-2)/n*ndiv)+1;
for j = 1:n-2,
	% Left side correlation
	C = xcorr1(prof2(profind(j),:), xsect(j,:));
	[maxcorr, ind] = max(C);
	
	lx(j+1) = mx(j+1) - dy(j)*off(ind);
	ly(j+1) = my(j+1) + dx(j)*off(ind);
	
	% Right side correlation (note the reversal of the profile)
	C = xcorr1(prof2(profind(j),end:-1:1), xsect(j,:));
	[maxcorr, ind] = max(C);
	
	rx(j+1) = mx(j+1) - dy(j)*off(ind);
	ry(j+1) = my(j+1) + dx(j)*off(ind);
end;

% Add on the edges of the image region so that all our variables correspond
% to the original image
mx = mx+rng(1)-1;
my = my+rng(3)-1;
lx = lx+rng(1)-1;
ly = ly+rng(3)-1;
rx = rx+rng(1)-1;
ry = ry+rng(3)-1;
