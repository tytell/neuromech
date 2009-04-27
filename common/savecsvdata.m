function savecsvdata(file,varargin)
% function savecsvdata(file,var1,var2,...)
%   or     savecsvdata(file,{var1,var2,...},{name1,name2,...})
%
% Options:
%   'varnames' - cell string containing the variable names
%   'nanstring' - String to output for a NaN value (default '.')
%   'emptystring' - String to output for an empty value (default '" "')
%   'flatten' - Flatten the variables before writing them to a file
%
% If no variable names are passed in, the function attempts to guess them
% using the names of the variables passed in as arguments

opt.varnames = {};
opt.nanstring = '.';
opt.emptystring = '" "';
opt.flatten = false;

[opt,data,inddata] = parsevarargin(opt,varargin,'firstoptionnumber',2,'leaveunknown');
nvar = length(data);

%handle the two old forms of parameters
if ((numel(data) == 1) && iscell(data{1}) && (numel(data{1}) > 1)),
    %savecsvdata(file,{var1,var2...})
    data = data{1};
elseif ((numel(data) == 2) && iscell(data{1}) && iscellstr(data{2}) && ...
        (length(data{1}) == length(data{2}))),
    %savecsvdata(file,{var1,var2,...},{name1,name2,...})
    opt.varnames = data{2};
    data = data{1};
end;

%try to get the names of the arguments if they didn't pass them in
if (isempty(opt.varnames) && (nargin >= length(data))),
    opt.varnames = cell(1,nvar);
    for i = 1:nvar,
        opt.varnames{i} = inputname(inddata(i)+1);
    end;
end;

if (opt.flatten),
    for i = 1:nvar,
        data{i} = data{i}(:);
    end;
else
    %break the data up into columns
    for i = 1:nvar,
        nc = size(data{i},2);
        
        data{i} = mat2cell(data{i},size(data{i},1),ones(1,nc));
        vn1 = opt.varnames{i};
        opt.varnames{i} = cell(1,nc);
        
        for j = 1:nc,
            opt.varnames{i}{j} = sprintf('%s%d',vn1,j);
        end;
    end;
    data = cat(2,data{:});
    opt.varnames = cat(2,opt.varnames{:});
end;

nvar = size(data,2);

%construct output templates
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
                      opt.varnames{i}));
    end;
    
    len = max(len,length(data{i}));
end;

%check for short data sets
for i = 1:nvar,
    if (length(data{i}) < len),
        warning(sprintf(['Variable %s is short.  Adding missing ' ...
                         'values.'],opt.varnames{i}));
        
        if (isnumeric(data{i})),
            data{i}(end+1:len) = NaN;
        elseif (ischar(data{i}))
            data{i}(end+1:len) = ' ';
        elseif (iscell(data{i})),
            [data{i}{end+1:len}] = deal(opt.emptystring);
        end;
    end;
end;

if (length(data) ~= length(opt.varnames)),
    error('Names and data should be the same size');
end;

fid = fopen(file,'w');
if (fid == -1),
    error('Couldn''t open file for output.');
end;

fprintf(fid,'%s',opt.varnames{1});
fprintf(fid,',%s',opt.varnames{2:end});
fprintf(fid,'\n');

for i = 1:len,
    for j = 1:nvar,
        if (isnumeric(data{j}(i))),
            if (isfinite(data{j}(i)) && isreal(data{j}(i))),
                fprintf(fid,tplt{j},data{j}(i));
            else
                fprintf(fid,'%s',opt.nanstring);
            end;
        elseif (iscell(data{j}(i))),
            if (~isempty(data{j}(i))),
                fprintf(fid,'%s',data{j}{i});
            else
                fprintf(fid,'%s',opt.emptystring);
            end;
        else
            if (data{j}(i) ~= ' '),
                fprintf(fid,tplt{j},data{j}(i));
            else
                fprintf(fid,'%s',opt.emptystring);
            end;
        end;
        
        if (j < nvar),
            fprintf(fid,',');
        end;
    end;
    fprintf(fid,'\n');
end;

fclose(fid);
