function fn = getfilenames(d)

k = find(d == '\');
if (~isempty(k)),
	path = d(1:k(end));
else
	path = [];
end;

files = dir(d);
if (length(files) > 0),
	for i = 1:length(files),
		fn{i,1} = [path files(i).name];
	end;
else
	fn = {};
end;
