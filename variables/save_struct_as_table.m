function save_struct_as_table(outname, S)

C = struct2cell(S);
colnames = fieldnames(S);

ncols = size(C,1);

C = flatten(C,2:ndims(C));
for i = 1:numel(C)
    C{i} = C{i}(:);
end
if size(C,2) > 1
    C2 = cell(ncols,1);
    for i = 1:ncols
        C2{i} = cat(1,C{i,:});
    end
    C = C2;
end

nrows = cellfun(@numel, C);
nrows = sum(nrows,2);

nd = cellfun(@ndims, C);
sz = zeros(max(nd),length(C));
for i = 1:nd
    sz(i,:) = cellfun(@(x) size(x,i), C);
end

if any(nrows ~= nrows(1)) || any(sz(:) ~= flatten(repmat(sz(:,1),[1 ncols])))
    error('All elements in the struct must have the same size');
end
nrows = nrows(1);

tplt = cell(1,ncols);
coltplt = cell(1,ncols);
D = cell(nrows,ncols);
for i = 1:ncols
    if ~iscell(C{i})
        D(:,i) = num2cell(C{i}(:));
    else
        D(:,i) = C{i}(:);
    end
    if ischar(C{i}) || iscell(C{i})
        tplt{i} = '%s';
    elseif isnumeric(C{i})
        tplt{i} = '%f';
    else
        error('Unrecognized type in struct');
    end
    coltplt{i} = '%s';
end

tplt = strjoin(tplt,',');
coltplt = strjoin(coltplt,',');

fid = fopen(outname, 'w');
fprintf(fid, [coltplt '\n'], colnames{:});

D = D';
fprintf(fid, [tplt '\n'], D{:});
fclose(fid);


