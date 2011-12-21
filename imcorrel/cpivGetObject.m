function varargout = cpivGetObject(objspec, fr)

if (isempty(objspec)),
	varargout = cell(1,nargout);
elseif (nargout == 1),
	if (iscell(objspec)),
		objspec = objspec{1};
	end;
	
	if (ischar(objspec)),
		if (strfind(lower(objspec),'avi')),
			obj = frame2im(aviread(objspec,fr));
		else
			obj = imread(objspec);
		end;
	elseif (isnumeric(objspec)),
		obj = objspec;
	else
		error('Unknown object specification');
	end;
	
	if (range(obj) > 1),
		varargout{1} = obj > 128;
	else
		varargout{1} = obj > 0.5;
	end;
elseif (nargout == 4),
	out = objspec(2:end);
	if (length(out) == 2),
		out{3} = zeros(size(out{1}));
		out{4} = zeros(size(out{1}));
	end;
	if (any(size(out{1}) == 1)),
		out{1} = shiftdim(out{1});
		out{2} = shiftdim(out{2});
		out{3} = shiftdim(out{3});
		out{4} = shiftdim(out{4});
	else
		out{1} = out{1}(:,fr);
		out{2} = out{2}(:,fr);
		out{3} = out{3}(:,fr);
		out{4} = out{4}(:,fr);
	end;
	varargout = out;
end;

