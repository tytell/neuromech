function varargout = resizeallvars(dim,good, varargin)
% function varargout = resizeallvars(dim,good, varargin)

opt.verbose = false;
opt.struct = struct([]);
[opt,vars] = parsevarargin(opt,varargin, 2, 'leaveunknown',true);

if (isempty(vars))
    if (~isempty(opt.struct))
        where = 'struct';
        S = opt.struct;
        vars = fieldnames(S);
        vals = struct2cell(S);
        sz = cellfun(@size,vals,'UniformOutput',false);
        sz = catuneven(1,sz{:});
    else
        vv = evalin('caller','whos');
        vars = {vv.name};
        sz = catuneven(1,vv.size);
        where = 'caller';
        vals = {};
    end;
elseif ((numel(vars) == 1) && ischar(vars{1}) && ...
        (strcmp(vars{1},'base') || strcmp(vars{1},'caller'))),
    where = vars{1};
    vv = evalin(vars{1},'whos');
    vars = {vv.name};
    sz = catuneven(1,vv.size);
    vals = {};
else
    where = 'caller';
    vals = vars;
    sz = cellfun(@size,vals, 'UniformOutput',false);
    sz = catuneven(1,sz{:});
end;
sz(isnan(sz)) = 1;
    
if (strcmp(where,'struct'))
    if (nargout ~= 1)
        error('Wrong number of outputs');
    else
        varargout = cell(1,1);
    end;
elseif (nargout ~= 0)
    if (nargout ~= length(vars))
        error('Wrong number of outputs');
    else
        varargout = cell(1,length(vars));
    end;
end;

canresize = sz(:,dim) == length(good);

ngood = sum(good);
for i = 1:length(vars),
    if (isempty(vals))
        V1 = evalin(where,vars{i});
    else
        V1 = vals{i};
    end;
    
    if (canresize(i))
        sz1 = sz(i,:);
        rest = [1:dim-1 dim+1:ndims(V1)];
        pmt = [dim rest];
        
        V1 = permute(V1,pmt);
        V1 = reshape(V1,[sz1(dim) prod(sz1(rest))]);
        V1 = V1(good,:);
        V1 = reshape(V1, [ngood sz1(rest)]);
        V1 = ipermute(V1,pmt);
        
        if (nargout == 0)
            assignin(where,vars{i},V1);
        elseif (strcmp(where,'struct'))
            S.(vars{i}) = V1;
        else
            varargout{i} = V1;
        end;
            
        if (opt.verbose)
            sz2 = size(V1);
            szstr1 = sprintf('%dx',sz1);
            szstr1 = szstr1(1:end-1);
            szstr2 = sprintf('%dx',sz2);
            szstr2 = szstr2(1:end-1);
            
            fprintf('%s: %s -> %s\n', vars{i}, szstr1, szstr2);
        end;
    elseif ((nargout ~= 0) && ~strcmp(where,'struct'))
        varargout{i} = vals{i};
    end;
end;

if (strcmp(where,'struct'))
    varargout{1} = S;
end;



