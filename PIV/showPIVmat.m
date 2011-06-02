function showPIVmat(files,sigmin,snrmin,stdvel)

if (nargin == 1),
	sigmin = 0.6;
	snrmin = 1;
	stdvel = 3;
end;

if (ischar(files)),
	files = {files};
end;

ff = get(0,'Children');
if (~isempty(ff)),
	B = ff(1);
else
	B = figure;
end;
if (length(ff) > 1),
	A = ff(2);
else
	A = figure;
end;

set(0,'CurrentFigure',A);
clf;
set(0,'CurrentFigure',B);
clf;

for i = 1:length(files),
	load(files{i},'x','y','u','v','DU','I1','data');
	
	dx = x(1,2) - x(1,1);
	dy = y(2,1) - y(1,1);
	

	mag = sqrt(u.^2 + v.^2);
	k = find(((data.Error ~= 0) & (data.Error ~= 16)) | (data.Signal < sigmin) | ...
			(data.SNR < snrmin) | (abs(mag-nanmean(mag(:))) > 3*nanstd(mag(:))));
	
	u(k) = NaN;
	v(k) = NaN;
	
	[us,vs,DU] = agw(x+u/2, y+v/2, u, v, x,y);
	w = DU.dudy - DU.dvdx;

	k = find((data.Error ~= 0) & (data.Error ~= 16));
	us(k) = NaN;
	vs(k) = NaN;
	w(k) = NaN;
	
	umed = nanmedian(u(:));
	
	set(0,'CurrentFigure',A);
	imshow(I1,'n');
	hold on;
	quiverc(x,y,u-umed,v,'y','truncate');
	set(gcf,'Renderer','painters');
	hold off;
	set(A,'Name',sprintf('%s unsmoothed',files{i}));
	
	set(0,'CurrentFigure',B);
	pcolor(x-dx/2,y-dy/2,w);
	shading interp;
	axis equal ij;
	hold on;
	quiverc(x,y,us-umed,vs,'k','l','truncate');
	colormap default; 
	
	c1 = min(w(:));
	c2 = max(w(:));
	cr = max(abs(c1),abs(c2));
	caxis([-cr cr]);
	colorbar;
	hold off;
	set(B,'Name',sprintf('%s smoothed',files{i}));
	
	if (i < length(files)),
		input('Press return');
	end;
end;
