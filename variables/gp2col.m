function [C,gpcol,m,n] = gp2col(x,varargin)
%GP2COL   Separates grouped data into columns
%   [C,gpcol,m,n] = gp2col(x,groups...)
%     or            gp2col(x,{gp1,gp2,...})
%   For matrix x and groups, takes each unique group (or combination of groups, 
%   if there are multiple groups) and turns it into a column in C.  
%   gpcol indicates which group that each column corresponds to.  It has
%   the same number of rows as groups, and the same number of columns as C.
%
%   m and n are indices that work similarly to the m and n output
%   parameters from UNIQUE.  Specifically, C(m) = x and x(n) = C (almost).
%   Since there may be different numbers of items in each group, there will
%   be values in m and n that are equal to zero.  So we have to write
%
%   good = m ~= 0;
%   C(m(good)) = x(good);
%
% See also UNIQUE and IM2COL.

% Mercurial revision hash: $Revision$ $Date$
% See https://bitbucket.org/tytell/matlab_variables/overview
% Copyright (c) 2011, Eric Tytell <tytell at jhu dot edu>

opt.quiet = false;
opt.keepnans = true;

%first parse the group matrix
%turn the group values (like 'A', 'B', 'C' or 5,7,9) into indices (1,2,3)
if (iscell(varargin{1})),
    gp = varargin{1};
    
    %cell array of things the same size as x
    gpind = zeros(numel(x), length(gp));
    gpvals = cell(length(gp),1);
    gpnames = cell(length(gp),1);
    for i = 1:length(gp),
        if ((ndims(gp) ~= ndims(x)) || any(size(gp) ~= size(x))),
            error('All group variables should be the same size as x.');
        end;
        [gpvals{i}, ~, gpind(isfinite(gp{i}),i)] = unique(gp{i}(isfinite(gp{i})));
        gpnames{i} = num2str(i);
    end;
    
    opts = 2;
else
    %check for multiple group variables passed directly
    isgp = false(size(varargin));
    for i = 1:length(varargin),
        if (((ndims(varargin{i}) == ndims(x)) && all(size(varargin{i}) == size(x))) || ...
            (any(size(x) == 1) && (numel(varargin{i}) == numel(x)))),
            isgp(i) = true;
        end;
    end;
    
    if (any(isgp)),
        gparg = find(isgp);
        ngp = length(gparg);
        
        gpnames = cell(ngp,1);
        gpvals = cell(ngp,1);
        gpind = zeros(numel(x),ngp);
        
        for i = 1:length(gparg),
            gp1 = varargin{gparg(i)};
            
            [gpvals{i}, ~, gpind(isfinite(gp1),i)] = unique(gp1(isfinite(gp1)));
            gpnames{i} = inputname(gparg(i) + 1);
            if (isempty(gpnames{i})),
                gpnames{i} = num2str(i);
            end;
        end;
        
        opts = gparg(end)+1;
    else
        gp1 = varargin{1};
        szgp = size(gp1);
        
        if ((ndims(gp1) == ndims(x)+1) && all(size(x) == szgp(1:end-1))),
            %multiple groups, with the number of groups in the last dimension
            gp1 = reshape(gp1,[numel(x) szgp(end)]);

            gpnames = cell(size(gp1,2),1);
            gpvals = cell(size(gp1,2),1);
            gpind = zeros(size(gp1));

            for i = 1:size(gp1,2),
                nonan = isfinite(gp1(:,i));
                [gpvals{i}, ~, gpind(nonan,i)] = unique(gp1(nonan,i));
                gpnames{i} = num2str(i);
            end;
        elseif (any(size(x) == 1) && (ndims(gp1) == 2) && any(szgp == length(x))),
            %multiple groups, with x a vector and gp1 a matrix.  Be nice about
            %vector transposition
            if (size(x,1) == 1),
                x = x';
            end;
            if (size(gp1,2) == length(x)),
                gp1 = gp1';
            end;

            gpnames = cell(size(gp1,2),1);
            gpvals = cell(size(gp1,2),1);
            gpind = zeros(size(gp1));

            for i = 1:size(gp1,2),
                nonan = isfinite(gp1(:,i));
                [gpvals{i}, ~, gpind(nonan,i)] = unique(gp1(nonan,i));
                gpnames{i} = num2str(i);
            end;
        else
            error('Group matrix must be the same size as x, or have its last dimension larger than x');
        end;
        
        opts = 2;
    end;
end;

opt = parsevarargin(opt,varargin(opts:end),opts);
        
ngp = size(gpind,2);
ngpvals = cellfun(@length,gpvals);

%now that we've sorted out the groups individually, find out how many
%unique combinations there are
good = all(gpind > 0,2);
N = numel(x);

combind = zeros(size(x));
[combvalind, ~, combind(good)] = unique(gpind(good,:), 'rows');
ncomb = size(combvalind,1);
ncombmax = prod(ngpvals);
if (any(~good) && opt.keepnans),
    ncomb = ncomb+1;
    combind(~good) = ncomb;
    combvalind(end+1,:) = NaN;
    isnangp = true;
else
    isnangp = false;
end;

if (~opt.quiet),
    fprintf('Found %d groups.\n', size(gpind,2));
    
    for i = 1:ngp,
        fprintf('Group %s values: ',gpnames{i});
        switch class(gpvals{i}),
            case 'logical',
                vals = char(size(gpvals{i}));
                vals(gpvals{i}) = 'T';
                vals(~gpvals{i}) = 'F';
                
                fprintf('%c ',vals);
                
            case 'char',
                fprintf('%c ',gpvals{i});
                
            case 'cell',
                if (iscellstr(gpvals{i})),
                    fprintf('%s ',gpvals{i}{:});
                elseif (all(cellfun(@isnumeric,gpvals{i}))),
                    fprintf('%g ',gpvals{i}{:});
                else
                    fprintf('Mixed type cell array with %d values.',length(gpvals{i}));
                end;
                
            otherwise,
                fprintf('%g ',gpvals{i});
        end;
        fprintf('\n');
    end;
    
    fprintf('\n%d combinations represented (%d%% of %d possible)\n', ...
        ncomb,round(ncomb/ncombmax*100),ncombmax);
    if (isnangp || ~opt.keepnans),
        fprintf('%d elements lost due to NaN group values\n', sum(~good));
    end;
end;
                
%now run through and find out how many times each combination is replicated 
ncombrep = zeros(1,ncomb);
row = zeros(N,1);
col = zeros(N,1);
for i = 1:N,
    if (combind(i) ~= 0),
        ncombrep(combind(i)) = ncombrep(combind(i)) + 1;
        row(i) = ncombrep(combind(i));
        col(i) = combind(i);
    end;
end;

%define the column matrix
C = NaN([max(ncombrep) ncomb]);
m = zeros(size(x));
m(good) = sub2ind(size(C),row(good),col(good));
C(m(good)) = x(good);

%and which combination each column corresponds to
if (all(cellfun(@isnumeric,gpvals))),
    gpcol = zeros(ngp,ncomb);
    b = ncomb;
    if (isnangp),
        b = b-1;
    end;
    for i = 1:b,
        for j = 1:ngp,
            gpcol(j,i) = gpvals{j}(combvalind(i,j));
        end;
    end;
    if (isnangp),
        gpcol(:,end) = NaN;
    end;
else
    gpcol = cell(ngp,ncomb);
    b = ncomb;
    if (isnangp),
        b = b-1;
    end;    
    for i = 1:b,
        for j = 1:ngp,
            gpcol{j,i} = gpvals{j}(combvalind(i,j));
        end;
    end;
    if (isnangp),
        [gpcol{:,end}] = deal(NaN);
    end;
end;

%and the index matrix to turn it back into x
xind = 1:numel(x);
n = zeros(size(C));
n(m(good)) = xind(good);



    

        
       