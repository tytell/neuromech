function saveSystat(file,data,names)

nvar = length(data);

len = 0;
for i = 1:nvar,
	if (ischar(data{i})),
		tplt{i} = '%10c';
	elseif (isnumeric(data{i})),
		tplt{i} = '%10g';
	else
		error(sprintf('Unsuported data type in variable %s.',names{i}));
	end;
	
	data{i} = data{i}(:);
	len = max(len,length(data{i}));
end;

for i = 1:nvar,
	if (length(data{i}) < len),
		warning(sprintf('Variable %s is short.  Adding missing values.',names{i}));
	
		if (isnumeric(data{i})),
			data{i}(end+1:len) = NaN;
		elseif (ischar(data{i}))
			data{i}(end+1:len) = ' ';
		end;
	end;
end;

fid = fopen(file,'w');
if (fid == -1),
	error('Couldn''t open file for output.');
end;

fprintf(fid,'%10s',names{1});
fprintf(fid,', %10s',names{2:end});
fprintf(fid,'\n');

for i = 1:len,
	for j = 1:nvar,
		if (isnumeric(data{j}(i))),
			if (isfinite(data{j}(i))),
				fprintf(fid,tplt{j},data{j}(i));
			else
				fprintf(fid,'%10s','.');
			end;
		else,
			if (data{j}(i) ~= ' '),
				fprintf(fid,tplt{j},data{j}(i));
			else
				fprintf(fid,'%10s','" "');
			end;
		end;
		
		if (j < nvar),
			fprintf(fid,', ');
		end;
	end;
	fprintf(fid,'\n');
end;

fclose(fid);
