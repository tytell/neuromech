function varargout = getvar(varargin)
%GETVAR   Loads variables out of the base workspace into a function
% getvar(variables...)
%   or     getvar(S,variables)
%   or     good = getvar(...)
%
%   Loads variables out of the base workspace into the calling function (or
%   from a file or structure)
%
% Options:
%    '-file',file - Loads variables from a file, much like 'load'
%    '-fromstruct',S - Loads variables from a struct S.  (also can use
%                   '-struct'
%    '-tostruct' - Returns a structure containing the requested variables.
%    '-all' - Operates on all of the current variables
%    '-except',variables - Except for these variables
%    '-keepundef' - Do not clear requested variables that weren't found
%
% Returns true if it found all of the requested variables.

% Mercurial revision hash: $Revision$ $Date$
% See https://bitbucket.org/tytell/matlab_variables/overview
% Copyright (c) 2011, Eric Tytell <tytell at jhu dot edu>

isvar = true(size(varargin));
exceptind = [];
optind = [];
infile = '';
usestruct = false;
ishidden = false;
isall = false;
clearundef = true;
warnundef = true;
tostruct = false;

%process options
if ((nargin >= 1) && isstruct(varargin{1})),
    %first parameter is a structure
    data = varargin{1};
    usestruct = true;
    isvar(1) = false;
    p = 2;
elseif ((nargin >= 1) && exist(varargin{1},'file')),
    %first parameter might be a file
    warning('getvar:oldsyntax',...
        'Possibly using old syntax to read from file.  Use -file instead');
    p = 1;
else
    p = 1;
end;

%look for options
while (p <= length(varargin)),
    if (~ischar(varargin{p})),
        error('getvar:badarg',...
            'Unrecognized argument at position %d.  Maybe quotes are missing?', p);
    end;
    switch lower(varargin{p}),
        case '-file',   
            infile = varargin{p+1};
            isvar(p:p+1) = false;
            optind(end+1) = p;
            p = p+2;
            
        case {'-struct','-fromstruct'},
            data = varargin{p+1};
            usestruct = true;
            isvar(p:p+1) = false;
            optind(end+1) = p;
            p = p+2;
            
        case '-tostruct',
            tostruct = true;
            isvar(p) = false;
            optind(end+1) = p;
            p = p+1;
            
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
            
        case '-keepundef',
            clearundef = false;
            optind(end+1) = p;
            isvar(p) = false;
            p = p+1;
            
        case '-nowarn',
            warnundef = true;
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
        k = length(varargin);
    end;
    
    except = varargin(exceptind+1:k);
    isvar(exceptind+1:k) = false;
else
    except = {};
end;

%pull out the requested variable names from varargin
requestedvars = varargin(isvar);

%find the existing variable names, based on the mode
if (usestruct),
    existingvars = fieldnames(data);
    errstr = ' in struct';
elseif (~isempty(infile) && exist(infile,'file') && ~exist(infile,'dir')),
    existingvars = who('-file',infile);
    errstr = sprintf(' in file %s',infile);
elseif (ishidden),
    if (~isappdata(0,'HiddenData')),
        warning('getvar:nohiddendata','There are no hidden variables');
        data = struct;
    else
        data = getappdata(0,'HiddenData');
    end;
    usestruct = true;
    errstr = ' in hidden data';
else
    existingvars = evalin('base','who');
    errstr = '';
end;

%handle -all or -except with no other args
if (isall || (isempty(requestedvars) && ~isempty(except))),
    requestedvars = existingvars;
end;

%remove exceptions
remove = ismember(requestedvars,except);
requestedvars = requestedvars(~remove);

%good requests are those that exist
good = ismember(requestedvars, existingvars);

if (~isempty(infile)),
    %we put this off, because we don't want to load the whole file if we
    %don't need to
    data = load(infile, requestedvars{good});
    usestruct = true;
end;

%save the requested variables that don't exist
undef = requestedvars(~good);
requestedvars = requestedvars(good);

if (usestruct),
    %run through the structure
    data = getfieldsonly(data, requestedvars);
    
    if (tostruct),
        outstruct = data;
    else
        for i = 1:length(requestedvars),
            assignin('caller',requestedvars{i},data.(requestedvars{i}));
        end;
    end;
else
    %or pull the variables out of the base workspace and assign them in the
    %caller
    if (~isempty(requestedvars))
        for i = 1:length(requestedvars),
            var = evalin('base',requestedvars{i});
            if (tostruct)
                outstruct.(requestedvars{i}) = var;
            else
                assignin('caller',requestedvars{i},var);
            end;
        end;
    else
        outstruct = struct;     % empty structure
    end;
end;

if (tostruct)
    if (nargout == 1),
        varargout = {outstruct};
    elseif (nargout == 2),
        varargout = {all(good), outstruct};
    end;
elseif (any(~good)),
    %if some variables were not defined, clear them in the caller, so that the
    %calling function doesn't use an old value
    if (clearundef),
        clearcmd = ['clear ' sprintf('%s ',undef{:})];
        clearcmd(end) = ';';
        
        evalin('caller',clearcmd);
    end;
    
    if (nargout == 0),
        if (warnundef),
            undef = sprintf('%s, ',undef{:});
            warning('getvar:undefinedVar','getvar requestedvars undefined variable(s) %s%s', ...
                undef(1:end-2), errstr);
        end;
    elseif (nargout == 1),
        varargout{1} = all(good);
    end;
elseif (nargout == 1),
    varargout{1} = true;
end;


