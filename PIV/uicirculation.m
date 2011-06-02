function varargout = uicirculation(x,y,u,v,varargin)
% function varargout = uicirculation(x,y,u,v,varargin)

type = 'circle';
def = [];
if ((nargin == 6) & all(size(varargin{1}) == size(varargin{2})) & isnumeric(varargin{1})),
	type = 'preset';
	xpt = shiftdim(varargin{1});
	ypt = shiftdim(varargin{2});
else,
	for i = 1:length(varargin),
		if (ischar(varargin{i})),
			type = varargin{i};
		else
			def = varargin{i};
		end;
	end;
end;

h = findobj('Tag','Circulation');
if (~isempty(h)),
	delete(h);
end;

dx = diff(x,[],2);
dy = diff(y,[],1);
dr = nanmean(abs([dx(:); dy(:)]));

if (strncmpi(type,'circle',length(type))),
	if (isempty(def)),
		[ctrx,ctry,maxr] = selectCircle;
	else
		ctrx = def(1);
		ctry = def(2);
		maxr = def(3);
	end;
	
	r = dr/2:dr:maxr;
	for i = 1:length(r),
		npt = ceil(2*pi*r(i)/dr);
		theta = linspace(0,2*pi,npt+1);
		xpt = ctrx + r(i).*cos(theta);
		ypt = ctry + r(i).*sin(theta);
		
		ur = interp2(x,y,u, xpt,ypt, '*cubic');
		vr = interp2(x,y,v, xpt,ypt, '*cubic');
		
		upar = -sin(theta).*ur + cos(theta).*vr;
		
		rcirc(i) = trapz(r(i)*theta, upar);
	end;
	
	sgn = median(sign(diff(rcirc)));
	rc = rcirc * sgn;
	rc = (rc - min(rc))/range(rc);
	
	if (isempty(def)),
		scale = max(abs(rcirc))/r(end);
		h = addplot(xpt,ypt,'k',ctrx+r, ctry+rcirc/scale, 'r-');
		set(h,'Tag','Circulation');
	end;
	
	if (nargout == 0),
		[c,ind] = last(rcirc);
		fprintf('Circulation = %f at radius %f\n',c,r(ind));
	else
		varargout = {rcirc,r};
	end;
elseif (strncmpi(type,'line',length(type))),
	if (isempty(def)),
		[lx,ly] = selectLine;
	else
		lx = def(1:2);
		ly = def(3:4);
	end;
	
	len = sqrt(diff(lx).^2 + diff(ly).^2);
	s = 0:dr:len;
	m = [diff(lx) diff(ly)]/len;
	
	xpt = lx(1) + m(1) * s;
	ypt = ly(1) + m(2) * s;
	
	upt = interp2(x,y,u, xpt,ypt, '*cubic');
	vpt = interp2(x,y,v, xpt,ypt, '*cubic');
	
	velpar = (diff(lx)*upt + diff(ly)*vpt)/len;
	nonan = find(isfinite(velpar));
	
	if (length(nonan) > 1),
		circ = trapz(s(nonan),velpar(nonan));
	else
		circ = NaN;
	end;
	scale = 0.2*len/max(abs(velpar));
	
	if (nargout == 0),
		h = addplot(xpt,ypt, 'k-', xpt - m(2)*velpar*scale, ypt + m(1)*velpar*scale, 'r-');
		set(h,'Tag','Circulation');
		
		fprintf('Circulation = %f along length %f\n', circ, len);
	end;
	
	varargout = {circ,xpt,ypt,velpar};
elseif (strncmpi(type,'preset',length(type))),
	npt = length(xpt);
	
	s = [0; cumsum(sqrt(diff(xpt).^2 + diff(ypt).^2))];
	m(2:npt-1,1) = (xpt(3:end)-xpt(1:end-2))./(s(3:end)-s(1:end-2));
	m(2:npt-1,2) = (ypt(3:end)-ypt(1:end-2))./(s(3:end)-s(1:end-2));
	m(1,1) = (xpt(2)-xpt(1))./s(2);
	m(1,2) = (ypt(2)-ypt(1))./s(2);
	m(npt,1) = (xpt(end)-xpt(end-1))./(s(end)-s(end-1));
	m(npt,2) = (ypt(end)-ypt(end-1))./(s(end)-s(end-1));

	upt = interp2(x,y,u, xpt,ypt, '*cubic');
	vpt = interp2(x,y,v, xpt,ypt, '*cubic');
	
	velpar = (m(:,1).*upt + m(:,2).*vpt);
	
	circ = trapz(s,velpar);
	varargout = {circ,xpt,ypt,velpar};
end;



	