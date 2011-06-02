function saveInsight(name, w,h, x,y,u,v,err)

fid = fopen(name,'w');
if (fid == -1)
	error(sprintf('Couldn''t open %s for writing.\n',name));
end;

fprintf(fid, 'TITLE = "%s" VARIABLES= "X pixel", "Y pixel", "U pixel", "V pixel", "CHC" ZONE T="Pixel, Height=%d, Width=%d " I=%d, J=%d, F=POINT\n', ...
					name, h,w,size(x,1),size(x,2));

x = x';
y = h - y';
u = u';
v = -v';
err = err';

err(err ~= 0) = -1;
err(err == 0) = 1;

fprintf(fid, '%f, %f, %f, %f, %d\n',[x(:) y(:) u(:) v(:) err(:)]');

fclose(fid);
