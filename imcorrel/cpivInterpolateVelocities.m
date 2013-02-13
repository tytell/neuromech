function [ug,vg] = cpivInterpolateVelocities(x,y,u,v,err, xp,yp, grid,fr, interpmode)

[bx,by, bu,bv] = cpivGetObject(grid.Object, fr);

switch (lower(grid.ObjInterp)),
case 'uv',
	[ug,vg] = doInterp(x,y,u,v,err, xp,yp, interpmode, {bx,by, bu,bv});
case 'nt',
	[enx,eny, q1,q2, bnx,bny] = normtang(x,y, bx,by);
	
	n = enx.*u + eny.*v;
	t = -eny.*u + enx.*v;
	bn = bnx.*bu + bny.*bv;
	bt = repmat(NaN, size(bn));		% don't restrict tangential velocities
	
	[ng,tg] = doInterp(x,y,n,t, err, xp,yp, interpmode, {bx,by, bn,bt});
	
	[enxp,enyp] = normtang(xp,yp, bx,by);
	
	ug = ng.*enxp - tg.*enyp;
	vg = ng.*enyp + tg.*enxp;
case 'none',
	[ug,vg] = doInterp(x,y,u,v, err, xp,yp, interpmode, {});
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% doInterp
function [ug,vg] = doInterp(x,y,u,v,err, xp,yp, mode, extra)

% tpaps smoothing value based on Spedding and Rignot
smooth = 0.0737;

switch lower(mode),
case 'agw',
	k = find(err ~= 0);
	u(k) = NaN;
	v(k) = NaN;
	[ug,vg] = agw(x,y, u,v, xp,yp, [], extra{:});
case 'sts',
	if (~isempty(extra)),
		x = [x(:); extra{1}(:)];
		y = [y(:); extra{2}(:)];
		u = [u(:); extra{3}(:)];
		v = [v(:); extra{4}(:)];
		err = [err(:); zeros(size(extra{1}(:)))];
	end;

	if (any(isnan(u) ~= isnan(v))),
		k1 = find(isfinite(u) & (err == 0));
		sp1 = tpaps([x(k1) y(k1)]', u(k1)', smooth);
		
		k2 = find(isfinite(v) & (err == 0));
		sp2 = tpaps([x(k2) y(k2)]', v(k2)', smooth);
		
		ug = repmat(NaN, size(xp));
		vg = repmat(NaN, size(yp));

		k = find(isfinite(xp) & isfinite(yp));
		uu = fnval(sp1, [xp(k) yp(k)]);
		ug(k) = uu;
		vv = fnval(sp2, [xp(k) yp(k)]);
		vg(k) = vv;
	else	
		k = find(isfinite(u) & isfinite(v) & (err == 0));
		
		sp = tpaps([x(k) y(k)]', [u(k) v(k)]', smooth);
		
		ug = repmat(NaN, size(xp));
		vg = repmat(NaN, size(yp));
		
		k = find(isfinite(xp) & isfinite(yp));
		uv = fnval(sp, [xp(k) yp(k)]');
		ug(k) = uv(1,:);
		vg(k) = uv(2,:);
	end;
otherwise,
	error(sprintf('Unrecognized interpolation mode ''%s''.', interpmode));
end;

