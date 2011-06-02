function [tab,combvals2,ncombvals,combind] = summarizegroups(varargin)
% function tab = summarizegroups(S)
%    or          summarizegroups(groups)
%    or          summarizegroups(gp1,gp2,gp3...)
%
%
% S is a structure, with each field representing a group.  Fields have to be the
% same size.  groups is a cell array, where each element in the array represents a group.
% Or groups can be passed as individual parameters (gp1, gp2, etc).
%
% Summarizes the number of elements in each combination of groups for multidimensional
% complex groupings.  Returns a table containing the number of elements, and also produces
% a (moderately) nice looking output.
%
% Options:
%   'groupnames' - Cell string array with the names of the groups.  If it's not passed in,
%   the groups are named for the field names (if a structure is passed in) or the names of
%   the input variables (if individual variables are passed in).

opt.gpnames = {};
opt.use = {};

if (isstruct(varargin{1})),
    S = varargin{1};
    opt.gpnames = fieldnames(S);
    groups = struct2cell(S);
    p = 2;
elseif (iscell(varargin{1})),
    groups = varargin{1};
    p = 2;
else
    p = 1;
    
    nd = ndims(varargin{1});
    sz = size(varargin{1});

    opt.gpnames = {};
    while ((p <= nargin) && (ndims(varargin{p}) == nd) && ...
            all(size(varargin{p}) == sz)),
        opt.gpnames{p} = inputname(p);
        p = p+1;
    end;
end;

opt = parsevarargin(opt,varargin(p:end),'firstoptionnumber',p);

if (~isempty(opt.use)),
    good = opt.use;
else
    good = true(size(groups{1}));
end;

ndall = cellfun(@ndims,groups);
nd = mode(ndall);
goodgp = ndall == nd;

szall = cellfun(@size,groups(goodgp),'UniformOutput',false);
szall = cat(1,szall{:});
[szu,q,szindall] = unique(szall,'rows');
szind = mode(szindall);
goodgp = szindall == szind;

if (any(~goodgp)),
    if (any(cellfun(@isempty,opt.gpnames(~goodgp)))),
        gpexc = sprintf('%d ',find(~goodgp));
    else
        gpexc = sprintf('%s ',opt.gpnames{~goodgp});
    end;
    
    warning('summarizegroups:differentgroupsizes',...
        'Some group variables are different sizes.  Excluding %s...', gpexc);
    
    groups = groups(goodgp);
    opt.gpnames = opt.gpnames(goodgp);
    szall = szall(goodgp,:);
end;

groups = cellfun(@(x) (x(:)), groups, 'UniformOutput',false);
good = good(:);

ngp = length(groups);
nvals = length(groups{1});

gpind = zeros(nvals,ngp);
gpvals = cell(1,ngp);
ngpvals = zeros(1,ngp);
for i = 1:ngp,
    if (isnumeric(groups{i})),
        nonan = ~isnan(groups{i});
        good1 = good & nonan;
    else
        nonan = true(size(groups{i}));
        good1 = good;
    end;
    
    [gpvals1,q,gpind1] = unique(groups{i}(good1));
    if (any(good & ~nonan)),
        gpvals1(end+1) = NaN;           %#ok
    end;
    
    gpvals{i} = gpvals1;
    gpind(good1,i) = gpind1;
    gpind(good & ~nonan,i) = length(gpvals1);

    ngpvals(i) = length(gpvals1);
end;

[combvals,q,combind] = unique(gpind,'rows');

ncomb = size(combvals,1);

combvals2 = cell(ncomb,ngp);
for i = 1:ngp,
    combvals2(:,i) = num2cell(gpvals{i}(combvals(:,i)));
end;

ncombvals = zeros(ncomb,1);
tab = zeros(ngpvals);
tabind = zeros(ngpvals);
for i = 1:ncomb,
    ncombvals(i) = sum(combind == i);
    
    c = num2cell(combvals(i,:));
    ind = sub2ind(ngpvals,c{:});
    
    tab(ind) = ncombvals(i);
    tabind(ind) = i;
end;

[ngpvals,ord] = sort(ngpvals,'descend');
ord([1 2]) = ord([2 1]);                % swap first two dimensions
ngpvals([1 2]) = ngpvals([2 1]);

gpvals = gpvals(ord);
opt.gpnames = opt.gpnames(ord);
combvals = combvals(:,ord);

tab2 = permute(tab,ord);
tab2ind = permute(tabind,ord);

oneval = find(ngpvals == 1);
if (~isempty(oneval)),
    fprintf('In all cases:\n    ');
    for i = 1:length(oneval),
        switch class(gpvals{oneval(i)}),
            case 'char',
                fprintf('%s=%c', opt.gpnames{oneval(i)},gpvals{oneval(i)});
            case 'cell',
                fprintf('%s=%s', opt.gpnames{oneval(i)},gpvals{oneval(i)}{:});
            otherwise,
                fprintf('%s=%g', opt.gpnames{oneval(i)},gpvals{oneval(i)});
        end;
        
        if (i < length(oneval)),
            fprintf(', ');
        else
            fprintf('\n');
        end;
    end;
end;

lastgood = last(ngpvals > 1);

tab2 = reshape(tab2,[ngpvals(1:2) prod(ngpvals(3:end))]);
tab2ind = reshape(tab2ind,size(tab2));
for i = 1:size(tab2,3),
    ind = tab2ind(:,:,i);
    if (any(ind(:) ~= 0)),
        ind = first(ind(:),ind(:) ~= 0);
        
        comb = combvals(ind,:);
        
        for j = 3:lastgood,
            switch class(gpvals{j}(comb(j))),
              case 'char',
                fprintf('%s=%c', opt.gpnames{j},gpvals{j}(comb(j)));
              case 'cell',
                fprintf('%s=%s', opt.gpnames{j},gpvals{j}{comb(j)});
              otherwise,
                fprintf('%s=%g', opt.gpnames{j},gpvals{j}(comb(j)));
            end;

            if (j < lastgood),
                fprintf(', ');
            else
                fprintf('\n');
            end;
        end;
        showtable(tab2(:,:,i),gpvals(1:2),opt.gpnames(1:2));
        fprintf('\n\n');
    end;
end;

combind = reshape(combind,szall(1,:));

function showtable(tab2d,gpvals2,gpnames2)

fprintf('%20s   %s ->\n',' ',gpnames2{2});
fprintf('%19s |  ',gpnames2{1});
switch class(gpvals2{2}(1)),
  case 'char',
    fprintf('%-10c ',gpvals2{2});
  case 'cell',
    if (ischar(gpvals2{2}{1})),
        fprintf('%-10s ',gpvals2{2}{:});
    end;
  otherwise,
    fprintf('%-10g ',gpvals2{2});
end;
fprintf('\n');
fprintf('%19s V\n',' ');

for i = 1:size(tab2d,1),
    switch class(gpvals2{1}(i)),
      case 'char',
        fprintf(' %20c: ',gpvals2{1}(i));
      case 'cell',
        if (ischar(gpvals2{1}{i})),
            fprintf(' %20s: ',gpvals2{1}{i});
        end;
      otherwise,
        fprintf(' %20g: ',gpvals2{1}(i));
    end;
    fprintf('%-10d ',tab2d(i,:));
    fprintf('\n');
end;
