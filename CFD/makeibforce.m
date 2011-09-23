function F = makeibforce(varargin)
% function IB = makeibforce(fxltot,fyltot,fxrtot,fyrtot,fxmtot,fymtot,fxntot,fyntot, ...
%    fxrmus,fyrmus,fxlmus,fylmus,fmusPr,fmusPl,actl,actr,...
%    fxsp,fysp,energy)
%   or          makeibforce('base')
%   or          makeibforce(S)
%
% Helper function to make an IBforce structure, which is required by some
% of the other functions.  Takes the coordinates of the forces on each of
% the immersed boundary lines, plus the forces due to the muscles and the
% activation parameters.  As long as the names of the parameters are the
% same as those above, the order doesn't matter.  Or it can read the
% parameters directly out of the base workspace ('base' option), or out of
% a structure (S).

Fnames = {'fxltot','fyltot','fxrtot','fyrtot','fxmtot','fymtot','fxntot','fyntot', ...
    'fxrmus','fyrmus','fxlmus','fylmus','fmusPr','fmusPl','actl','actr',...
    'fxsp','fysp','energy'};

n = nargin;
if ((n >= 2) && islogical(varargin{n}) && ~strcmp(inputname(n),'actl') && ...
        ~strcmp(inputname(n),'actr'))
    good = varargin{n};
    n = n-1;
end;

if ((n == 1) && ischar(varargin{1}) && ...
        strcmp(varargin{1},'base')),
    [ok,F] = getvar(Fnames{:},'-tostruct');
    if (~ok),
        error('Standard IB variables are not in base workspace');
    end;
elseif ((n == 1) && isstruct(varargin{1}))
    [ok,F] = getvar('-struct',varargin{1},Fnames{:},'-tostruct');
    if (~ok),
        warning('makeibforce:missingvars','Some IB variables are not in the structure');
    end;
else
    nm = cell(1,n);
    ok = true;
    for i = 1:n,
        nm{i} = inputname(i);
        ok = ~isempty(nm{i});
    end;
    
    if (ok && all(ismember(nm,Fnames))),
        F = cell2struct(varargin(1:n),nm,2);
    else
        F = cell2struct(varargin(1:n),Fnames(1:n),2);
    end;
end;

sz1 = structfun(@(x) (size(x,1)), F);
sz2 = structfun(@(x) (size(x,2)), F);

if (exist('good','var'))
    if ((size(good,1) == 1) && all(sz2 == size(good,2))),
        F = structfun(@(x) (x(:,good)), F, 'UniformOutput',false);
    elseif ((size(good,2) == 1) && all(sz1 == size(good,1)))
        F = structfun(@(x) (x(good,:)), F, 'UniformOutput',false);
    elseif (all(sz1 == size(good,1)) && all(sz2 == size(good,2)))
        F = structfun(@(x) (x(good)), F, 'UniformOutput',false);
    end;
end;
    
    