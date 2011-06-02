function varargout = putvar(varargin)
% function putvar(variables...)
%   or     S = putvar(S,variables)
%
% Options:
%    '-file',file - Saves variables to a file, like the save function,
%    except that it always appends to the file
%    '-struct',S - Puts variables into the structure S, not the base
%    workspace.  Allows functions based around putvar and getvar to be
%    converted into batch processing algorithms
%    '-all' - Operates on all of the current variables
%    '-except',variables - Except for these variables

vars = varargin;
isvar = true(size(vars));
exceptind = [];
optind = [];
outfile = '';
usestruct = false;
ishidden = false;
isall = false;

%process options
p = 1;
while (p <= length(vars)),
    switch lower(vars{p}),
        case '-file',   
            outfile = vars{p+1};
            isvar(p:p+1) = false;
            optind(end+1) = p;
            p = p+2;
            
        case '-struct',
            if ((p+1 <= length(vars)) && isstruct(vars{p+1}))
                data = vars{p+1};
                isvar(p:p+1) = false;     
                p = p+2;
            else
                data = struct;
                isvar(p) = false;
                p = p+1;
            end;
            usestruct = true;
            optind(end+1) = p;
            
        case '-all',
            isall = true;
            isvar(p) = false;
            optind(end+1) = p;
            p = p+1;
            
        case '-except',
            exceptind = p;
            isvar(p) = false;
            p = p+1;
            
        case '-hidden',
            ishidden = true;
            optind(end+1) = p;
            isvar(p) = false;
            p = p+1;
            
        otherwise,
            isvar(p) = true;
            p = p+1;
    end;
end;

%find exceptions
if (~isempty(exceptind)),
    if (any(optind > exceptind)),
        k = min(optind(optind > exceptind))-1;
    else
        k = length(vars);
    end;
    
    except = vars(exceptind+1:k);
    isvar(exceptind+1:k) = false;
else
    except = {};
end;

%process all variables, if necessary
if (isall),
    vars = evalin('caller','who');
else
    vars = vars(isvar);
end;

%remove excepted variables
remove = ismember(vars,except);
vars = vars(~remove);
    
if (usestruct),
    %save to a structure
    for i = 1:length(vars),
        data.(vars{i}) = evalin('caller',vars{i});
    end;
    
    if (nargout == 1),
        varargout = {data};
    end;
elseif (~isempty(outfile)),
    %save to a file
    for i = 1:length(vars),
        F.(vars{i}) = evalin('caller',vars{i});
    end;
    
    save(outfile,'-append','-struct','F');
elseif (ishidden),
    %save to hidden data
    if (isappdata(0,'HiddenData')),
        data = getappdata(0,'HiddenData');
    else
        data = struct;
    end;
    
    for i = 1:length(vars),
        var = evalin('caller',vars{i});
        data.(vars{i}) = var;
    end;
    setappdata(0,'HiddenData',data);
else
    %save to the base workspace
    for i = 1:length(vars),
        var = evalin('caller',vars{i});
        assignin('base',vars{i},var);
    end;
end;
