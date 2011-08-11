function varargout = putvar(varargin)
%PUTVAR   Puts variables from a function into the base workspace
% putvar(variables...)
%   or     putvar(S,fields)
%   or     S = putvar(...)
%
% Options:
%    '-file',file - Saves variables to a file, like the save function,
%    except that it always appends to the file
%    '-tostruct',S - Puts variables into the structure S, not the base
%    workspace.  Allows functions based around putvar and getvar to be
%    converted into batch processing algorithms
%    '-fromstruct',S - Puts variables in the structure S into the base
%    workspace.  (equivalent to putvar(S,fields...))
%    '-all' - Operates on all of the current variables
%    '-except',variables - Except for these variables
%
% Note that the syntax is a bit different than SAVE.  For example,
%   putvar('-fromstruct',S,'-all')
% saves all of the fields in the struct S to the base workspace, but
%   save('file.mat','-struct','S')
% (with the 'S' in quotes) saves the fields of S to a file.
%
% SEE ALSO
%   GETVAR, SAVE

% Mercurial revision hash: $Revision$ $Date$
% See https://bitbucket.org/tytell/matlab/overview
% Copyright (c) 2011, Eric Tytell <tytell at jhu dot edu>

vars = varargin;
isvar = true(size(vars));
exceptind = [];
optind = [];
outfile = '';
tostruct = false;
ishidden = false;
isall = false;
fromstruct = false;

p = 1;
if ((nargin >= 1) && isstruct(varargin{1})),
    %first parameter is a structure
    indata = varargin{1};
    fromstruct = true;
    isvar(1) = false;
    p = 2;
end;

%process options
while (p <= length(vars)),
    switch lower(vars{p}),
        case '-file',   
            outfile = vars{p+1};
            isvar(p:p+1) = false;
            optind(end+1) = p;
            p = p+2;
            
        case '-tostruct',
            if ((p+1 <= length(vars)) && isstruct(vars{p+1}))
                outdata = vars{p+1};
                isvar(p:p+1) = false;   
                outdataname = '';
                p = p+2;
            elseif ((p+1 <= length(vars)) && ischar(vars{p+1}) && (vars{p+1}(1) ~= '-'))
                outdata = struct;
                outdataname = vars{p+1};
                isvar(p:p+1) = false;
                p = p+2;
            else
                outdata = struct;
                outdataname = '';
                isvar(p) = false;
                p = p+1;
            end;
            tostruct = true;
            optind(end+1) = p;
           
        case '-fromstruct',
            indata = varargin{p+1};
            isvar(p:p+1) = false;
            p = p+2;
            fromstruct = true;
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
    if (fromstruct)
        vars = fieldnames(indata);
    else
        vars = evalin('caller','who');
    end;
else
    vars = vars(isvar);
end;


%remove excepted variables
remove = ismember(vars,except);
vars = vars(~remove);

if (ishidden)
    if (isappdata(0,'HiddenData')),
        outdata = getappdata(0,'HiddenData');
    else
        outdata = struct;
    end;
end;

for i = 1:length(vars)
    if (fromstruct)
        val = indata.(vars{i});
    else
        val = evalin('caller',vars{i});
    end;
    
    if (tostruct),
        %save to a structure
        outdata.(vars{i}) = val;
    elseif (~isempty(outfile)),
        %save to a file
        F.(vars{i}) = val;
    elseif (ishidden),
        outdata.(vars{i}) = val;
    else
        %save to the base workspace
        assignin('base',vars{i},val);
    end;
end;

if (tostruct && ~isempty(outdataname))
    assignin('caller',outdataname,outdata);
end;

if (tostruct),
    %save to a structure
    if (nargout == 1),
        varargout = {outdata};
    end;
elseif (~isempty(outfile)),
    save(outfile,'-append','-struct','F');
elseif (ishidden),
    setappdata(0,'HiddenData',outdata);
end;

