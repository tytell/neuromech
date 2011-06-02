function [enx,eny,perps,perpdist, bnx,bny] = normtang(x,y,bx,by)
% function [enx,eny,perps,perpdist, bnx,bny] = normtang(x,y,bx,by)
% Returns normal vectors (enx,eny) to the curve (bx,by) at the positions 
% (x,y).  Also returns the distance along the curve perps and the normal
% distance to the curve perpdist for each point (x,y).  bnx and bny
% are the normal vectors to the curve at each curve point (bx,by).
%
% Can take multiple curves, in which case different columns of (bx,by)
% refer to differnt curves.  Then (enx,eny) both are of size [size(x)
% size(bx,2)].

bx = shiftdim(bx);
by = shiftdim(by);

% set up the sizes right
enx = repmat(NaN,[size(x) size(bx,2)]);
eny = repmat(NaN,[size(y) size(by,2)]);
perpdist = repmat(NaN,[size(x) size(bx,2)]);
perps = repmat(NaN,[size(x) size(bx,2)]);
nxy = prod(size(x));

% remove NaNs
nonan = find(isfinite(x) & isfinite(y));
x = x(nonan)';
y = y(nonan)';

npt = length(nonan);

% check for closed curves
if ((bx(1) == bx(end)) & (by(1) == by(end))),
	isclosed = 1;
else
	isclosed = 0;
end;

nb = size(bx,1)-1;

% get dx/ds and dy/ds along the curve
dy = diff(by);
dx = diff(bx);
len = sqrt(diff(bx).^2 + diff(by).^2);
dy = dy./len;
dx = dx./len;

% calculate normal vectors at the curve itself
s = [zeros(1,size(bx,2)); cumsum(len)];
bnx = -deriv(s,by);
bny = deriv(s,bx);

timedWaitBar(0, 'Calculating normals and tangents.');

% iterate through all the curves
for i = 1:size(bx,2),
	bx1 = bx(:,i);
	by1 = by(:,i);
	dx1 = dx(:,i);
	dy1 = dy(:,i);
	len1 = len(:,i);
	
	% calculate distance between each (x,y) and each (bx1,by1)
	dist2 = (repmat(bx1(1:end-1),[1 npt]) - repmat(x, [nb 1])).^2 + ...
			(repmat(by1(1:end-1),[1 npt]) - repmat(y, [nb 1])).^2;
	
	% find the boundary point closest to each spatial point
	[md,ind] = min(dist2);

	% calculate the perpendicular distance to the line segment connecting 
	% the closest point and the next point along the curve.  dx1=0 or
	% dy1=0 are degenerate cases and have to be treated separately
	k = find((dx1(ind) ~= 0) & (dy1(ind) ~= 0));
	perpdist1(k,1) = (1./dx1(ind(k)).*(x(k)' - bx1(ind(k))) - 1./dy1(ind(k)).*(y(k)' - by1(ind(k)))) ./ ...
				(dx1(ind(k))./dy1(ind(k)) + dy1(ind(k))./dx1(ind(k)));
	perps1(k,1) = 1./dx1(ind(k)).*(x(k)' - bx1(ind(k))) - perpdist1(k).*dy1(ind(k))./dx1(ind(k));

	% treat degenerate cases
	x0 = find(dx1(ind) == 0);
	perpdist1(x0) = len1(ind(x0)).*(x(x0)' - bx1(ind(x0))) ./ dy1(ind(x0));
	y0 = find(dy1(ind) == 0);
	perpdist1(y0) = -len1(ind(y0)).*(y(y0)' - by1(ind(y0))) ./ dx1(ind(y0));
	
	% perps1 is now the intersection of the normal from each (x,y) to
	% the line segment between the closest point and the next one.  If
	% perps1 < 0, we should look at the previous line segment.  If it's
	% greater than the length of the line segment, we should look at the
	% next one
	k1 = find(perps1 < 0);
	k2 = find(perps1 > len1(ind));
	
	% adjust indices appropriately
	ind(k1) = ind(k1)-1;	
	ind(k2) = ind(k2)+1;

	% deal with wrapping around the end of the curve, if the curve is
	% closed
	if (isclosed)
		ind(ind == 0) = nb;
		ind(ind > nb) = 1;
	end;
	
	% make sure we don't try to look off the end of the curve
	q = find((ind < 1) | (ind > nb));
	ind(q) = 1;		% dummy value
	
	% and recalculate perpdist1 and perps1
	k = find((dx1(ind) ~= 0) & (dy1(ind) ~= 0));
	perpdist1(k) = (1./dx1(ind(k)).*(x(k)' - bx1(ind(k))) - 1./dy1(ind(k)).*(y(k)' - by1(ind(k)))) ./ ...
				(dx1(ind(k))./dy1(ind(k)) + dy1(ind(k))./dx1(ind(k)));
	perps1(k) = 1./dx1(ind(k)).*(x(k)' - bx1(ind(k))) - perpdist1(k).*dy1(ind(k))./dx1(ind(k));

	x0 = find(dx1(ind) == 0);
	perpdist1(x0) = len1(ind(x0)).*(x(x0)' - bx1(ind(x0))) ./ dy1(ind(x0));
	y0 = find(dy1(ind) == 0);
	perpdist1(y0) = -len1(ind(y0)).*(y(y0)' - by1(ind(y0))) ./ dx1(ind(y0));
	
	% points beyond the curve get NaNs
	perpdist1(q) = NaN;
	perps1(q) = NaN;
	
	% set the overall variable accorting to the nonan.  Make sure perps
	% refers to *total* distance along the curve, not just distance along
	% a line segment
	perps(nonan + (i-1)*nxy) = s(ind,i) + perps1;
	perpdist(nonan + (i-1)*nxy) = perpdist1;
	
	enx(nonan + (i-1)*nxy) = -dy1(ind) .* sign(perpdist1);
	eny(nonan + (i-1)*nxy) = dx1(ind) .* sign(perpdist1);
	
	if (~timedWaitBar(i/size(bx,2))),
		break;
	end;
end;

timedWaitBar(1);