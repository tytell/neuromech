function [x,y,h] = ginputseg(n)
% function [x,y,h] = ginputseg(n)
% Allows the user to select a line consisting of n points (n-1 segments).
% Returns the x and y coordinates of the line, plus the handles of the
% line segments.  If n is not passed, it allows unlimited segments.
% The user can stop before n points by pressing return, or cancel all points
% by pressing Esc.  They can also use the zoom function by pressing z, or zoom
% all the way out by pressing f.

% Set up n and the x,y, and h arrays.  It is more memory efficient to set up
% the number of points in the array at the beginning, if we know n.
if (nargin == 0)
	n = Inf;
elseif (isfinite(n))
	x = zeros(n,1);
	y = zeros(n,1);
	h = zeros(n,1);
end;

if (n == Inf)
	x = zeros(20,1);
	y = zeros(20,1);
	h = zeros(20,1);
end;

i = 1;
while (i <= n),
	% get the point
	[x1,y1,butt] = ginput(1);
	
	% check the button they pressed
	if (isempty(butt))		% return
		break;
	elseif (butt == 1)		% left mouse
		x(i,1) = x1;
		y(i,1) = y1;
		beep;
		if (i > 1)			% if we already have points, draw a line
			h(i) = line(x((i-1):i),y((i-1):i));
		else
			hold on;		% otherwise draw a point
			h(i) = plot(x1,y1,'+'); 
			hold off;
		end;
	elseif ((butt == 2) | (butt == 8))		% delete or right mouse
		if (i >= 0)
			delete(h(i-1));		% erase whatever we drew last
		end;
		i = i-2;				% and back up (it's -2 because we add 1 below)
	elseif (butt == 27)				% escape
		x = [];				% clear everything
		y = [];
		if (i > 1)
			delete(h(1:(i-1)));
		end;
		h = [];
		return;				% and exit
	elseif ((butt == 'c') | (butt == 'C'))		% close the curve
		if (i > 1),
			x(i) = x(1);
			y(i) = y(1);
			h(i) = line(x((i-1):i),y((i-1):i));
			break;
		end;
	elseif ((butt == 'z') | (butt == 'Z'))		% zoom
		i = i-1;	% we subtract one so that we won't record this point
		zoom on;
		input('Press return to continue');
	elseif ((butt == 'f') | (butt == 'F'))		% zoom out
		i = i-1;
		zoom out;
	end;
	
	i = i+1;
end;

% if we terminated before n points, make the size of x,y, and h right
i = i-1;
if (i < n)
	x = x(1:i);
	y = y(1:i);
	h = h(1:i);
end;

if (nargout == 1)
	x = [x y];
end;
