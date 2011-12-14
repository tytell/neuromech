function IB = makeib(varargin)
% function IB = makeib(xl,yl,xr,yr,xm,ym,xn,yn,ul,vl,ur,vr,um,vm,un,vn)
%   or          makeib('base')
%   or          makeib(S)
%
% Helper function to make an IB structure, which is required by some of the
% other functions.  Takes the coordinates of the l,r,m, and n lines and
% their velocities and produces a structure with the right components.  As
% long as the names of the parameters are the same as those above, the
% order doesn't matter.  Or it can read the parameters directly out of the
% base workspace ('base' option), or out of a structure (S).

IBnames = {'xl','yl','xr','yr','xm','ym','xn','yn',...
    'ul','vl','ur','vr','um','vm','un','vn'};

n = nargin;

if (n == 0),
    disp('No arguments to makeib.  If you want the base workspace, use the ''base'' option');
    return;
end;

if ((n >= 2) && islogical(varargin{n}))
    good = varargin{n};
    n = n-1;
end;

if ((n == 1) && ischar(varargin{1}) && ...
        strcmp(varargin{1},'base')),
    [ok,IB] = getvar(IBnames{:},'-tostruct');
    if (~ok),
        error('Standard IB variables are not in base workspace');
    end;
elseif ((n == 1) && isstruct(varargin{1}))
    [ok,IB] = getvar('-struct',varargin{1},IBnames{:},'-tostruct');
    if (~ok),
        warning('makeib:missingvars','Some IB variables are not in the structure');
    end;
else
    nm = cell(1,n);
    ok = true;
    for i = 1:n,
        nm{i} = inputname(i);
        ok = ~isempty(nm{i});
    end;
    
    if (ok && all(ismember(nm,IBnames))),
        IB = cell2struct(varargin,nm,2);
    else
        IB = cell2struct(varargin,IBnames(1:n),2);
        if (n < length(IBnames))
            warning('makeib:missingvars','Some IB variables do not appear to be in the arguments');
        end;
    end;
end;

sz1 = structfun(@(x) (size(x,1)), IB);
sz2 = structfun(@(x) (size(x,2)), IB);

if (exist('good','var'))
    if ((size(good,1) == 1) && all(sz2 == size(good,2))),
        IB = structfun(@(x) (x(:,good)), IB, 'UniformOutput',false);
    elseif ((size(good,2) == 1) && all(sz1 == size(good,1)))
        IB = structfun(@(x) (x(good,:)), IB, 'UniformOutput',false);
    elseif (all(sz1 == size(good,1)) && all(sz2 == size(good,2)))
        IB = structfun(@(x) (x(good)), IB, 'UniformOutput',false);
    end;
end;
    
    