function [x,y, wsz,ssz] = cpivGetPoints(gd, fr, pass)

wsz = gd.WindowSize(pass,:);
ssz = gd.SearchSize(pass,:);

[x,y] = meshgrid(gd.Start(pass,1) + ceil(wsz(1)/2) + ...
					(0:(gd.NPt(pass,1)-1)) * gd.Offset(pass,1), ...
				gd.Start(pass,2) + ceil(wsz(2)/2) + ...
					(0:(gd.NPt(pass,2)-1)) * gd.Offset(pass,2));

% deal with the object, if there is one
if (~isempty(gd.Object)),
	obj = cpivGetObject(gd.Object, fr);
	
	[wndx,wndy] = meshgrid((1:wsz(1)) - floor(wsz(1)/2), ...
							(1:wsz(2)) - floor(wsz(2)/2));
	
	nwnd = length(wndx(:));
	npt = length(x(:));
	
	wndx = repmat(x(:)', [nwnd 1]) + repmat(wndx(:), [1 npt]);
	wndy = repmat(y(:)', [nwnd 1]) + repmat(wndy(:), [1 npt]);
	
	wndobj = obj(sub2ind(size(obj), wndy, wndx));
	nobj = sum(wndobj)/nwnd;
	
	x(nobj > gd.Overlap) = NaN;
	y(nobj > gd.Overlap) = NaN;
end;

	
	