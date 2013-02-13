function [x,y] = cpivCheckGrid(width,height, sWidth,sHeight, x,y)

sw = ceil(sWidth/2);
sh = ceil(sHeight/2);

[i,j] = find(any(x-sw < 1, 3) | any(x+sw > width, 3) | ...
				any(y-sh < 1, 3) | any(y+sh > height, 3));

ng = size(x,3);
if (~isempty(i)),
	ind = sub2ind(size(x), repmat(i,[1 ng]), repmat(j,[1 ng]), ...
					repmat(1:ng, [length(i) 1]));
	
	x(ind) = NaN;
	y(ind) = NaN;
end;