function inside = findpivshadow(x,y, bx,by, dir, extra)
% function inside = findpivshadow(x,y, bx,by, dir, extra) 
%
% Finds all the points in x,y that should be shadowed by a boundary defined
% by bx,by.  bx,by should be mostly horizontal or vertical.  Convex regions
% will produce unpredictable results.  By default, it assumes that the
% shadow region is where the x,y values are lower than those in bx,by.
% In principle, dir = -1 changes this behavior, but I've never tested it.
%
% x,y - Meshgrid style array of points.  Size [m n]
% bx,by - Array of boundaries.  Size [npt nfr].
% dir - Either 1 (default) or -1.
%
% Returns an index value for all the points that should be shadowed,
% assuming that there is a variable somewhere with size [m n nfr].
%
% Example:
%
% Given x,y (size [68 78]), u,v (size [68 78 100]), and ox,oy (size [20
% 100])
%
% k = findpivshadow(x,y, ox,oy);
% u(k) = NaN;  v(k) = NaN;

if (numel(dir) ~= 2),
    error('dir must be a two element vector, specifying the shadow direction');
end;

%first rotate so that dir is to the right (x positive)
dir = dir ./ sqrt(sum(dir).^2);     % make dir a unit vector
xrot = dir(1)*x - dir(2)*y;
yrot = dir(2)*x + dir(1)*y;
bxrot = dir(1)*bx - dir(2)*by;
byrot = dir(2)*bx + dir(1)*by;

xmax = max(xrot(:));

inside = false(size(x,1),size(x,2),size(bx,2));
for i = 1:size(bx,2),
    bx1 = bxrot(:,i);
    by1 = byrot(:,i);
    
    bx1 = [xmax; bx1; xmax; xmax];
    by1 = by1([1 1:end end 1]);
    
    inside(:,:,i) = inpolygon(xrot,yrot, bx1,by1);
end;
