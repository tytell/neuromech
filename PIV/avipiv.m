function avipiv(aviname, varargin)
% function avipiv(aviname, out, grid, frames, ...)
% Runs superpiv on frames from an avi.  Writes the data to matlab files specified
% by out.  out can either be a cell string array of each output file name,
% corresponding to each frame, or a single string, to which the frame number of
% the first image will be attached.  grid is the PIV grid, as normally defined
% for superpiv.  frames is a list of frame numbers, either a vector containing the
% first frame number for each image pair (the second image will be the next frame),
% or a two column matrix, in which the first column contains the first image's frame
% number and the second column contains its pair's frame number.
%
% You can also pass any options that you can pass to superpiv.
%
% In the matlab output files, it saves the following variables:
%		x,y,u,v - The vector field calculated by superpiv
%		DU - The derivative matrix from superpiv
%		data - The data structure from superpiv
%		frame1, frame2, aviname - The frame numbers operated on in aviname
%		I1, I2 - The actual images
%		grid, options - The grid and options passed to superpiv

k = 1;
if (ischar(varargin{k}) & strfind(varargin{k},'.avi')),
	avibname = varargin{k};
	k = k+1;
else
	avibname = '';
end;
if (ischar(varargin{k})),
	out = varargin{k};
	k = k+1;
else
	out = [aviname(1:end-4) '.mat'];
end;
if (~isnumeric(varargin{k}) | ...
		((size(varargin{k},2) ~= 3) & (size(varargin{k},2) ~= 6)))
	error('You must pass a grid argument');
end;
grid = varargin{k};
k = k+1;

frames = shiftdim(varargin{k});
k = k+1;

options = varargin(k:end);

info = aviinfo(aviname);
avilen = info.NumFrames;

if ((~iscell(frames) & ~isnumeric(frames)) | (iscell(frames) & ~isnumeric(frames{1}))),
	error('You must pass an argument containing the frames to do PIV on.');
elseif (iscell(frames)),
	df = frames{1};
	frames = (df:df:avilen)';
elseif (size(frames,2) == 1),
	if ((min(frames) < 1) | (max(frames)+1 > avilen)),
		error(sprintf('Vector frame numbers must be between 1 and %d.',avilen-1));
	end;
else
	if ((min(frame(:)) < 1) | (max(frames(:)) > avilen)),
		error(sprintf('Matrix frame numbers must be between 1 and %d.',avilen));
	end;
end;

F = [];
k = strmatch('single',options);
if (~isempty(k)),
	isManyOutFiles = 0;
	options = options([1:k-1 k+1:end]);
	
	fid = fopen(out,'r');
	if (fid ~= -1),				% file exists
		fclose(fid);
		F = load(out);
		if (any(size(F.grid) ~= size(grid)) | any(F.grid ~= grid)),
			error('Output file exists with different grid.  Please delete or rename it.');
		end;
		
		k = [1 find(aviname == '\')+1];
		avif1 = aviname(k(end):end);
		k = [1 find(F.aviname == '\')+1];
		avif2 = aviname(k(end):end);
		
		if (~strcmpi(avif1,avif2)),
			error('Output file exists for different AVI.  Please delete or rename the output file.');
		end;
	end;
else
	isManyOutFiles = 1;
end;

overwriteOption = 0;
j = 1;
for i = 1:size(frames,1),
	frame1 = frames(i,1);
	if (size(frames,2) == 2)
		frame2 = frames(i,2);
	elseif (~isempty(avibname)),
		frame2 = frame1;
	elseif (i < avilen),
		frame2 = frame1+1;
	else
		break;
	end;
	overwrite = 1;
	
	if (iscellstr(out))
		outname = out{i};
	else,
		outname = sprintf('%s%04d.mat',out,frame1);
	end;

	if (isManyOutFiles),
		fid = fopen(outname,'r');
		if (fid ~= -1),
			overwrite = [];
			fclose(fid);
		end;
	end;
	
	if (~isempty(F) & (any(F.frames(:,1) == frame1))),
		overwrite = [];
	end;
	
	if (isempty(overwrite)),
		if (overwriteOption == 0),
			disp('Frame already calculated in existing output file.');
			ans = input('Overwrite? (y, n, a = always, b = never) ','s');
			
			switch (lower(ans)),
			case 'y',
				overwrite = 1;
			case 'n',
				overwrite = 0;
			case 'a',
				overwrite = 1;
				overwriteOption = 1;
			case 'b',
				overwrite = 0;
				overwriteOption = 2;
			otherwise,
				warning('Unrecognized option.  Not overwriting.');
				overwrite = 0;
			end;
		elseif (overwriteOption == 1),
			overwrite = 1;
		elseif (overwriteOption == 2),
			overwrite = 0;
		end;
	end;
	
	if (overwrite),
		if (isempty(avibname)),
			mov = aviread(aviname, [frame1 frame2]);
			I1 = frame2im(mov(1));
			I2 = frame2im(mov(2));
		else
			I1 = frame2im(aviread(aviname,frame1));
			I2 = frame2im(aviread(avibname,frame2));
		end;
		
		[x,y,u,v,DU,data] = superpiv(I1,I2,grid,'file',[i-1 size(frames,1)],options{:});
	
		if (isManyOutFiles),
			save(outname,'x','y','u','v','DU','data','aviname','avibname','frame1','frame2','I1','I2', ...
						'grid','options','-mat');
		else
			if (~isempty(F)),
				q = find(F.frames(:,1) == frame1);
				F.frames(q,:) = NaN;
			end;
			
			um(:,:,j) = u;
			vm(:,:,j) = v;
			DUm(j,1) = DU;
			datam(j,1) = data;
			framesm(j,:) = [frame1 frame2];
			
			j = j+1;
		end;
	end;
end;

% superpiv('end');

if (~isManyOutFiles),
	if (~isempty(F)),
		k = find(isfinite(F.frames(:,1)));
		u = cat(3,F.u(:,:,k),um);
		v = cat(3,F.v(:,:,k),vm);
		DU = [shiftdim(F.DU(k)); DUm];
		data = [shiftdim(F.data(k)); datam];
		frames = [F.frames(k,:); framesm];
	else
		u = um;
		v = vm;
		DU = DUm;
		data = datam;
		frames = framesm;
	end;
	
	[ff,sortorder] = sort(frames(:,1));
	u = u(:,:,sortorder);
	v = v(:,:,sortorder);
	DU = DU(sortorder);
	data = data(sortorder);
	frames = frames(sortorder,:);
	
	save(out,'x','y','u','v','grid','options','frames','data','DU','aviname','avibname','-mat');
end;


