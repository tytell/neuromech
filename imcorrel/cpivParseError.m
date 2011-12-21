function cpivParseError(error)

error = sort(error(:));
k = find(diff(error) > 0);
error = error([k; end]);
error = error(error > 0);

titles = {'Error Value:  ', 'L','R','T','B','Corr','Hart'};
type = {'d','s','s','s','s','s','s'};
len = num2cell(cellfun('size',titles,2) + 2);

maketplt = {len{:}; type{:}};
tplt = sprintf('%%-%d%s',maketplt{:});
tplt = [tplt '\n'];

titletplt = sprintf('%%-%ds',len{:});
fprintf([titletplt '\n'], titles{:});
fprintf([repmat('-',[1 sum(vertcat(len{:}))]) '\n']);

xes = cell(0);
for i = 1:length(error),
	bits = bitget(error(i),1:6);
	on = find(bits == 1);
	off = find(bits == 0);
	
	xes(on) = {'x'};
	xes(off) = {' '};
	
	fprintf(tplt,error(i),xes{:});
end;



