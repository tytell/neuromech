function [u,v] = tpapspiv(phi, smooth)

[dx,dy] = meshgrid((1:size(phi,2))-(size(phi,2)+1)/2, (1:size(phi,1))-(size(phi,1)+1)/2);

for i = 1:size(phi,3),
	for j = 1:size(phi,4),
		phi1 = phi(:,:,i,j);
		if (all(phi1 == phi1(1,1)) | any(isnan(phi1))),
			u(i,j) = NaN;
			v(i,j) = NaN;
		else
			sp = tpaps([dx(:) dy(:)]',phi1(:)',smooth);
			
			[uv,val,flag] = fminsearch(@spfun, [0 0], [], sp, [dx(1) dx(end) dy(1) dy(end)]);
			if (flag > 0),
				u(i,j) = uv(1);
				v(i,j) = uv(2);
			else
				u(i,j) = NaN;
				v(i,j) = NaN;
			end;
		end;
	end;
end;

function f = spfun(x, sp, rgn)
if ((x(1) < rgn(1)) | (x(1) > rgn(2)) | (x(2) < rgn(3)) | (x(2) > rgn(4))),
	f = 0;
else
	f = -fnval(sp, x);
end;
