function [gd, params] = cpivProcessGridParams(params, w,h, rgn)

grid = params{1};
params = params(2:end);

if (size(grid,2) == 3),
	grid = grid(:,[1 1 2 2 3 3]);
end;

% reset width and height according to our region
if (~isempty(rgn)),
	w = rgn(2)-rgn(1);
	h = rgn(4)-rgn(3);
end;

% minimum corner position
ax = ceil(grid(:,3)/2 - grid(:,1)/2) + 1;
ay = ceil(grid(:,4)/2 - grid(:,2)/2) + 1;

% number of grid squares (with fractional part)
nx = (w-1 - 2*ax-1 - (grid(:,1) - grid(:,5)))./grid(:,5);
ny = (h-1 - 2*ay-1 - (grid(:,2) - grid(:,6)))./grid(:,6);

% shift starting points according to the region
if (~isempty(rgn)),
	ax = ax+rgn(1)-1;
	ay = ay+rgn(3)-1;
end;

% and add on half the excess to center the grid
ax = ax + round(grid(:,5).*(nx - floor(nx))/2);
ay = ay + round(grid(:,6).*(ny - floor(ny))/2);

% remove fractional part of nx and ny
nx = floor(nx);
ny = floor(ny);

gd.Start = [ax ay];
gd.NPt = [nx ny];
gd.WindowSize = grid(:,1:2);
gd.SearchSize = grid(:,3:4);
gd.Offset = grid(:,5:6);
gd.Region = rgn;

gd.NPasses = size(grid,1);


