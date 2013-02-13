function [velperp,velpar,s,ut,vt, xt,yt] = transect(x,y,u,v, scale, x1,y1)
% function [velperp,velpar,s,ut,vt, xt,yt] = transect(x,y,u,v, scale, x1,y1)

if (nargin < 7),
	x1 = [];
	y1 = [];
	if (nargin < 5),
		scale = [];
	end;
end;

if (isempty(x1)),
	h = findobj(gca,'Tag','Transect');
	if (~isempty(h)),
		delete(h);
	end;
	
	[x1,y1] = selectline;
	
	interact = 1;
else
	interact = 0;
end;

dx = diff(x,[],2);
dy = diff(y,[],1);
d = nanmedian([dx(:); dy(:)]);

len = sqrt(diff(x1).^2 + diff(y1).^2);
m = [diff(x1); diff(y1)]/len;
n = floor(len/d);

s = 0:d:d*n;
xt = x1(1) + m(1)*s;
yt = y1(1) + m(2)*s;

ut = interp2(x,y,u, xt,yt);
vt = interp2(x,y,v, xt,yt);

velperp = (diff(y1)*ut - diff(x1)*vt)/len;
velpar = (diff(x1)*ut + diff(y1)*vt)/len;

if (interact),
	if (isempty(scale)),
		mv = max([abs(velperp) abs(velpar)]);
		ax = axis;
		if (abs(m(1)) > abs(m(2))),
			d = ax(4)-ax(3);
		else
			d = ax(2)-ax(1);
		end;
		scale = 0.1*d/mv;
	end;

	h = addplot(xt,yt,'k:',xt+velperp*m(2)*scale,yt-velperp*m(1)*scale,'r-',...
				xt-velpar*m(2)*scale,yt+velpar*m(1)*scale,'k--', ...
				'LineWidth',2);
	set(h,'Tag','Transect');
end;
