function calibGrid(I,dist,ang)
% function [d,dd,n, dx,ddx,nx, dy,ddy,ny] = calibGrid(I,dist,[ang1 ang2])
% Calculates the grid spacing based on an image of a
% rectangular grid.  The grid can be at an angle.
% Uses a Hough (radon) transform to determine the
% angle of the grid, then locates the lines and estimates
% their spacing.  The optional argument dist gives the distance
% between grid lines to calculate a scale factor.  The optional 
% argument ang specifies a range of angles to test, from a minimum 
% to a maximum, or 0 to do no rotation.
% 
% When spacing has been calculated, it returns the spacing
% as the median value, with the error as 0.7413 times the
% interquartile range (non-biased estimate for sigma based
% on a normal distribution).  We use median and IQR to avoid
% the effect of outliers.

% tolerence for integer spacings
tol = 0.1;

if (nargin < 3),
	ang = NaN;
	if (nargin < 2),
		dist = NaN;
	end;
end;

% First normalize the image so that it ranges from 0 to 1
I = im2double(I);
I = 1 - I;      % then invert it

% Display the image
imshow(I);
sz = size(I);

% Select the grid area
disp('Zoom to the grid');
zoom on;
input('Press return to continue');
ax = round(axis);

% Restrict the image to the selected area
I = I(ax(3):ax(4), ax(1):ax(2));

%%%% Calculate the grid angle

% Initial, coarse step, angles from -45 to 45 in steps
% of 5 deg.  Or if ang was passed, in 20 steps from
% the minimum to the maximum in ang.

% Radon is the Matlab implementation of a linear Hough
% transform.  Peaks in the radon transform indicate the
% presence of lines.
if (ang == 0),
	fprintf(1,'No rotation\n');
else
	if (isnan(ang))
		theta = -45:5:45;
	else
		theta = linspace(min(ang),max(ang),20);
	end;

	[R,pos] = radon(I, theta);
	
	% We should have a series of peaks at a single angle.
	% The sum of the radon transform across the rows is by
	% definition equal to the sum of all the elements in I,
	% so we exponentiate to pick out the peak.
	peak = sum(exp(R));
	[m,i] = max(peak);
	ang = theta(i);
	
	% Second, fine step.  Look +-2.5 deg from ang in steps
	% of 0.05 deg.
	theta = (ang-2.5):0.05:(ang+2.5);
	[R,pos] = radon(I, theta);
	
	peak = sum(exp(R));
	[m,i] = max(peak);
	ang = theta(i);
	
	fprintf(1,'Angle appears to be %g\n',ang);
	
	% Now rotate the image so that the grid is not at an angle
	I = imrotate(I, -ang, 'bilinear');
end;
imshow(I);

x = 1:size(I,2);
y = 1:size(I,1);

% The peaks corresponding to the grid lines will show up
% if we sum the image across the x and y directions
xpeak = sum(I);
ypeak = sum(I');

% To find the peaks, we find the zeros of the (appr) derivative
% of xpeak
[xz,s] = findzero(1.5:size(I,2), diff(xpeak));
xz = xz(find(s > 0));        % we only want maxima
hold on;
plot(x,xpeak/max(xpeak)*size(I,1)/2,'r:');
for i = 1:length(xz),
	hx(i) = plot([xz(i) xz(i)],[y(1); y(end)],'r--');
end;
hold off;

% Same as for x
[yz,s] = findzero(1.5:size(I,1), diff(ypeak));
yz = yz(find(s > 0));
hold on;
plot(ypeak/max(ypeak)*size(I,2)/2,y,'g:');
for i = 1:length(yz),
	hy(i) = plot([x(1) x(end)],[yz(i); yz(i)],'g--');
end;
hold off;

kx = 1:length(xz);
ky = 1:length(yz);

disp('Click to remove elements');
[rmx,rmy] = ginput(1);
while (~isempty(rmx)),
	[mx,xi] = min(abs(rmx - xz(kx)));
	[my,yi] = min(abs(rmy - yz(ky)));
	
	if (mx < my),
		delete(hx(kx(xi)));
		kx = kx([1:xi-1 xi+1:end]);
	else
		delete(hy(ky(yi)));
		ky = ky([1:yi-1 yi+1:end]);
	end;
	
	[rmx,rmy] = ginput(1);
end;

% Rough grid spacing is the median difference between the peak positions
rawdx = diff(xz(kx));
dx = median(rawdx);

% Some lines may be multiple grid lines apart, so correct for that
nspacex = rawdx/dx;
kx = find(abs(mod(nspacex+0.5,1) - 0.5) < tol);		% find integer spacings
nspacex = round(nspacex);
alldx = rawdx(kx)./nspacex(kx);
dx = mean(alldx);
nx = length(kx);                     % count
ddx = std(alldx)/sqrt(nx);

scalex = dx/dist;
dscalex = ddx/dist;

rawdy = diff(yz(ky));
dy = median(rawdy);
nspacey = rawdy/dy;
ky = find(abs(mod(nspacey+0.5,1) - 0.5) < tol);		% find integer spacings
nspacey = round(nspacey);
alldy = rawdy(ky)./nspacey(ky);
dy = mean(alldy);
ny = length(ky);                     % count
ddy = std(alldy)/sqrt(ny);

scaley = dy/dist;
dscaley = ddy/dist;

% Overall mean
d = mean([alldx alldy]);
n = length(kx) + length(ky);
dd = std([alldx alldy])/sqrt(n);

scale = d/dist;
dscale = dd/dist;

fprintf('dx = %g+-%g pix, n = %g\ndy = %g+-%g pix, n = %g\n', dx,ddx,nx, dy,ddy,ny);
fprintf('Overall: d = %g+-%g pix, n = %g\n', d,dd,n);

if (isfinite(scale)),
	fprintf('xscale = %g+-%g, yscale = %g+-%g\n',scalex,dscalex, scaley,dscaley);
	fprintf('scale = %g+-%g\n', scale,dscale);
end;

% Test whether the x and y means are different
Sp = sqrt(((nx-1)*std(alldx)^2 + (ny-1)*std(alldy))/(nx+ny-2));
t = (dx - dy)/(Sp*sqrt(1/nx + 1/ny));
P = 1 - tcdf(t,nx+ny-2);
if (P < 0.05)
	fprintf(1,'Warning: x and y means are different (t = %g, P = %g)\n',t,P);
else
	fprintf(1,'x and y means are the same (t = %g, P = %g)\n',t,P);
end;

