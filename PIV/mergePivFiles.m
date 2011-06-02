function mergePivFiles(files, out)

if (ischar(files)),
	files = getfilenames(files);
end;

F = load(files{1});
x = F.x;
y = F.y;
grid = F.grid;
options = F.options;
aviname = F.aviname;

u(:,:,1) = F.u;
v(:,:,1) = F.v;
frames(1,:) = [F.frame1 F.frame2];
data(1) = F.data;
DU(1) = F.DU;

for i = 2:length(files),
	F = load(files{i});
	if (~isfield(F,'x') | ~isfield(F,'y') | ~isfield(F,'u') | ~isfield(F,'v')),
		warning(sprintf('File %s is not a vector file.  Skipping it...', files{i}));
	elseif (~isfield(F,'aviname') | ~strcmpi(aviname,F.aviname)),
		warning(sprintf('File %s was from a different avi or not from an avi.  Skipping it...', files{i}));
	elseif ((ndims(F.x) ~= ndims(x)) | (any(size(F.x) ~= size(x))) | any(F.x ~= x) | ...
			(ndims(F.y) ~= ndims(y)) | (any(size(F.y) ~= size(y))) | any(F.y ~= y))
		warning(sprintf('File %s has a different grid.  Skipping it...', files{i}));
	else
		u(:,:,i) = F.u;
		v(:,:,i) = F.v;
		frames(i,:) = [F.frame1 F.frame2];
		data(i) = F.data;
		DU(i) = F.DU;
	end;
end;

[ff sortorder] = sort(frames(:,1));
frames = frames(sortorder,:);
u = u(:,:,sortorder);
v = v(:,:,sortorder);
data = data(sortorder);
DU = DU(sortorder);

save(out,'x','y','u','v','grid','options','frames','data','DU','aviname','-mat');

