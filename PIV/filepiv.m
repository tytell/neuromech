function filepiv(files, out, grid, varargin)

options = varargin;
N = size(files,1);

I1 = imread(files{1});

obj = [];
bx = [];
by = [];
bu = [];
bv = [];

k = [];
i = 1;
while (i <= length(options)),
	if ((ndims(options{i}) == 3) & all(size(options{i}) == [size(I1) N])),
		obj = options{i};
		k = [k i];
	elseif ((i < length(options)) & isnumeric(options{i}) & ...
			any(size(options{i}) == N) & ...
			all(size(options{i}) == size(options{i+1}))),
		bx = options{i};
		by = options{i+1};
		k = [k i i+1];
		if ((i+3 <= length(options)) & ...
				all(size(options{i+2}) == size(options{i})) & ...
				all(size(options{i+3}) == size(options{i}))),
			bu = options{i+2};
			bv = options{i+3};
			k = [k i+2 i+3];
			i = i+2;
		end;
		i = i+1;
	end;
	i = i+1;
end;
options = options(setdiff(1:length(options),k));

for i = 1:N,
	if (size(files,2) == 2),
		file1 = files{i,1};
		file2 = files{i,2};
	elseif (i < N),
		file1 = files{i};
		file2 = files{i+1};
	else
		break;
	end;
	
	I1 = imread(file1);
	I2 = imread(file2);

	objopt = {};
	if (~isempty(obj)),
		objopt = {objopt{:} obj(:,:,i)};
	end;
	if (~isempty(bx)),
		k = find(isfinite(bx(:,i)) & isfinite(by(:,i)));
		objopt = {objopt{:} bx(k,i) by(k,i)};
		if (~isempty(bu)),
			objopt = {objopt{:} bu(k,i) bv(k,i)};
		end;
	end;
	[x,y,u,v,DU,data] = superpiv(I1,I2,grid,'file',[i-1 length(files)],options{:},objopt{:});
	
	save(out{i},'x','y','u','v','DU','data','file1','file2','grid','options',...
				'I1','I2','-mat');
end;
