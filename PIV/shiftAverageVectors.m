function [x0,y0, ushift,vshift,wshift] = shiftAverageVectors(ctrx,ctry, ...
                                                  x,y,u,v,w)
% function [x0,y0, ushift,vshift,wshift] = 
%                       shiftAverageVectors(ctrx,ctry, x,y,u,v,w)

if (nargin == 6),
    w = [];
    isW = false;
else
    isW = true;
end;

ctrx0 = mean(ctrx);
ctry0 = mean(ctry);

if (all(x(1,:,1) == x(2,:,1))),
	isplaid = 1;
else
	isplaid = 0;
end;

if (isplaid),
	dx = x(1,2) - x(1);
	dy = y(2,1) - y(1);
	rx = [x(1) x(end)];
	ry = [y(1) y(end)];
else
	dx = diff(x(:,:,1),[],2);
	dx = nanmedian(dx(:));
	dy = diff(y(:,:,1),[],1);
	dy = nanmedian(dy(:));
	
	rx = [nanmedian(min(x(:,:,1),[],2)) nanmedian(max(x(:,:,1),[],2))];
	ry = [nanmedian(min(y(:,:,1),[],1)) nanmedian(max(y(:,:,1),[],1))];
end;

rx0 = rx - ctrx0;
ry0 = ry - ctry0;
rx0(1) = rx0(1) - mod(rx0(1),dx);
ry0(1) = ry0(1) - mod(ry0(1),dy);

[x0,y0] = meshgrid(rx0(1):dx:rx0(2), ry0(1):dy:ry0(2));

if (isplaid),
    if (isW),
        k = find(isnan(u) | isnan(v) | isnan(w));
    else
        k = find(isnan(u) | isnan(v));
    end;

    u(k) = 0;
    v(k) = 0;
    if (isW),
        w(k) = 0;
    end;

    %save the NaN indices so we can set them back to NaN after we shift
    %everything
    [r0,c0,fr0] = ind2sub(size(u),k);
end;

tr = round((ctry-y(1))/dy) + 1;
tc = round((ctrx-x(1))/dx) + 1;
tr0 = round(-ry0(1)/dy) + 1;
tc0 = round(-rx0(1)/dx) + 1;

ushift = repmat(NaN,[size(x0) length(ctrx)]);
vshift = repmat(NaN,[size(y0) length(ctry)]);
if (isW),
    wshift = repmat(NaN,[size(x0) length(ctrx)]);
else
    wshift = [];
end;

rx = sort(rx);
ry = sort(ry);

timedWaitBar(0,'Interpolating...');
for fr = 1:length(ctrx),
	if (size(x,3) > 1),
		x1 = x(:,:,fr) - ctrx(fr);
		y1 = y(:,:,fr) - ctry(fr);
	else
		x1 = x - ctrx(fr);
		y1 = y - ctry(fr);
	end;
	
	i = find((y0(:,1) >= ry(1)-ctry(fr)) & (y0(:,1) <= ry(2)-ctry(fr)));
	j = find((x0(1,:) >= rx(1)-ctrx(fr)) & (x0(1,:) <= rx(2)-ctrx(fr)));
	
	if (isplaid),
		u1 = repmat(NaN,size(x0));
		v1 = repmat(NaN,size(y0));
		u1(i,j) = interp2(x1,y1,u(:,:,fr), x0(i,j),y0(i,j), '*spline');
		v1(i,j) = interp2(x1,y1,v(:,:,fr), x0(i,j),y0(i,j), '*spline');
		if (isW),
                    w1 = repmat(NaN,size(x0));
                    w1(i,j) = interp2(x1,y1,w(:,:,fr), x0(i,j),y0(i,j), ...
                                      '*spline');
                end;

		q = find(fr0 == fr);
		r1 = r0(q) - tr(fr) + tr0;
		c1 = c0(q) - tc(fr) + tc0;
		
		q = find((r1 >= 1) & (r1 <= size(x0,1)) & ...
                         (c1 >= 1) & (c1 <= size(x0,2)));
		if (~isempty(q)),
                    k = sub2ind(size(u1),r1(q),c1(q));

                    u1(k) = NaN;
                    v1(k) = NaN;
                    if (isW),
                        w1(k) = NaN;
                    end;
		end;
		
		ushift(i,j,fr) = u1(i,j);
		vshift(i,j,fr) = v1(i,j);
                if (isW),
                    wshift(i,j,fr) = w1(i,j);
                end;
	else
		[u1,v1] = agw(x1,y1, u(:,:,fr),v(:,:,fr), x0(i,j),y0(i,j));
		ushift(i,j,fr) = u1;
		vshift(i,j,fr) = v1;
	end;
	
	timedWaitBar(fr/length(ctrx));
end;
