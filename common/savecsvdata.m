function saveSystat(file,data,names)

nvar = length(data);

sz(1,:) = cellfun('size',data,1);
sz(2,:) = cellfun('size',data,2);

if (all(sz(1,:) == sz(1,1)) & any(sz(2,:) ~= sz(2,1))),
    len = sz(1,1);

    k = 1;
    for i = 1:nvar,
        for j = 1:sz(2,i),
            d{k+j-1} = data{i}(:,j);
        end;
        k = k+sz(2,i);
    end;
    nvar = size(d,2);
    data = d;
end;

len = 0;
for i = 1:nvar,
    if (ischar(data{i})),
        tplt{i} = '%c';
    elseif (isnumeric(data{i})),
        tplt{i} = '%g';
    elseif (iscellstr(data{i})),
        tplt{i} = '%s';
    else
        error(sprintf('Unsuported data type in variable %s.', ...
                      names{i}));
    end;
    
    data{i} = data{i}(:);
    len = max(len,length(data{i}));
end;

for i = 1:nvar,
    if (length(data{i}) < len),
        warning(sprintf(['Variable %s is short.  Adding missing ' ...
                         'values.'],names{i}));
        
        if (isnumeric(data{i})),
            data{i}(end+1:len) = NaN;
        elseif (ischar(data{i}))
            data{i}(end+1:len) = ' ';
        end;
    end;
end;

if (length(data) ~= length(names)),
    error('Names and data should be the same size');
end;

fid = fopen(file,'w');
if (fid == -1),
    error('Couldn''t open file for output.');
end;

fprintf(fid,'%s',names{1});
fprintf(fid,',%s',names{2:end});
fprintf(fid,'\n');

for i = 1:len,
    for j = 1:nvar,
        if (isnumeric(data{j}(i))),
            if (isfinite(data{j}(i)) & isreal(data{j}(i))),
                fprintf(fid,tplt{j},data{j}(i));
            else
                fprintf(fid,'%s','.');
            end;
        elseif (iscell(data{j}(i))),
            fprintf(fid,'%s',data{j}{i});
        else,
            if (data{j}(i) ~= ' '),
                fprintf(fid,tplt{j},data{j}(i));
            else
                fprintf(fid,'%s','" "');
            end;
        end;
        
        if (j < nvar),
            fprintf(fid,',');
        end;
    end;
    fprintf(fid,'\n');
end;

fclose(fid);
