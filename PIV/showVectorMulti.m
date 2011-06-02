function [montage, mov] = showVectorMulti(x,y, u,v, outx,outy, back, vperc,cperc,ratio, mode)

if (nargin < 11),
	mode = 1;
	if (nargin < 10),
		ratio = 1.5;
		if (nargin < 9),
			cperc = 0.5;
			if (nargin < 8),
				vperc = 2;
			end;
		end;
	end;
end;

nfr = size(u,3);
nc = ceil(sqrt(ratio*nfr));
nr = ceil(nc/ratio);

cscale = max(prctile(back(:),100-cperc),-prctile(back(:),cperc));
d = x(1,2) - x(1);
velscale = 0.9*d/prctile(sqrt(u(:).^2 + v(:).^2),100-vperc);
veltrunc = 100-0.75*vperc;

colormap default;
cmap = colormap;
cmap(1,:) = [1 1 1];
colormap(cmap);

set(gcf,'Color','w');

if (mode == 1),
	for i = 1:nfr,
		imagesc(x(1,:),y(:,1),back(:,:,i));
		caxis([-cscale cscale]);
		ax = axis;
		hold on;
		quiverc(x,y, u(:,:,i),v(:,:,i), 'k', 'noheads', 'l', -velscale, 't',veltrunc);
		if (size(outx,1) == 1),
			plot(outx(:,i),outy(:,i),'ko','LineWidth',2);
		elseif (~isempty(outx)),
			plot(outx(:,i),outy(:,i),'k-','LineWidth',2);
		end;
		hold off;
		axis(ax);
		
		set(gca, 'XTick', [], 'YTick', []);
		
		mov(i) = getframe;
	end;
	
	I1 = im2double(frame2im(mov(1)));
	iim = 1:size(I1,1);
	jim = 1:size(I1,2);
	
	montage = zeros(size(I1,1)*nr, size(I1,2)*nc, 3);
	
	k = 1;
	for i = 1:nr,
		for j = 1:nc,
			if (k <= nfr),
				montage(iim + (i-1)*length(iim), jim + (j-1)*length(jim),:) = im2double(frame2im(mov(k)));
			end;
			k = k+1;
		end;
	end;
	
	imshow(montage);
else
	w = 0.8/nc;
	h = 0.8/nr;
	yax = 0.9 - h;
	
	clf;
	
	k = 1;
	for i = 1:nr,
		xax = 0.1;
		for j = 1:nc,
			if (k <= nfr),
				axes('Position',[xax yax w h]);
				
				imagesc(x(1,:),y(:,1),back(:,:,i));
				caxis([-cscale cscale]);
				ax = axis;
				hold on;
				quiverc(x,y, u(:,:,k),v(:,:,k), 'k', 'noheads', 'l', -velscale, 't',veltrunc);
				if (size(outx,1) == 1),
					plot(outx(:,k),outy(:,k),'ko','LineWidth',2);
				elseif (~isempty(outx)),
					plot(outx(:,k),outy(:,k),'k-','LineWidth',2);
				end;
				hold off;
				axis equal ij off;
				axis(ax);
			end;
			xax = xax + w;
		end;
		yax = yax - h;
		k = k + 1;
	end;
end;

