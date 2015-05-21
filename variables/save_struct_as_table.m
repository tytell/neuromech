function save_struct_as_table(outname, S)

C = struct2cell(S);
colnames = fieldnames(S);

nrows = cellfun(@numel, C);
ncols = length(C);

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


