function [opts,param,ismatch] = matchQuivercOption(opts,param,name,nick,def)
% Mercurial revision hash: $Revision: dd099ded8eec $ $Date: 2010/08/11 13:58:29 $
% Copyright (c) 2010, Eric Tytell

if (nargin == 4),
    isdef = 0;
else
    isdef = 1;
end;
if (isdef && isnumeric(def)),
    needsNum = 1;
else
    needsNum = 0;
end;

q = find(cellfun('isclass',param, 'char'));
ind = strmatch(lower(name),lower(param(q)),'exact');
if (isempty(ind)),
    ind = strmatch(lower(nick),lower(param(q)),'exact');
end;
ind = q(ind);

if (isdef),
    if (isempty(ind)),
        if (~isfield(opts,name)),
            opts.(name) = def;
        end;
        ismatch = 0;
    else
        if (ind == length(param)),
            error('Option %s requires a parameter.',name);
        elseif (needsNum && ~isnumeric(param{ind+1})),
            error('Option %s requires a numeric parameter.',name);
        end;
        opts.(name) = param{ind+1};
        ismatch = 1;
        
        param = param([1:ind-1 ind+2:end]);
    end;
else
    if (isempty(ind)),
        if (~isfield(opts,name)),
            opts.(name) = 0;
        end;
        ismatch = 0;
    else
        opts.(name) = 1;
        param = param([1:ind-1 ind+1:end]);
    end;
end;
