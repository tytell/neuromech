function varargout = use(S,varargin)
% function [ok,missingvars] = use(S,varnames,options)
% 
% Puts variables from the structure S into the caller's workspace.
% varnames optionally lists which variables.  If it's not passed, then
% we'll use all of the fields in S.  Normally use clears the variables in 
% the caller's workspace before it assigns them.  That way, if it fails on
% one, an old version won't exist.
%
% Options:
%     '-noclear' - Don't clear variables in the caller's workspace.
%
% Output: ok - boolean for whether it found all the variables.
%         missingvars - cellstr listing the variables it didn't find.
%
% Credit: Idea from JLab by J.M. Lilly, but code is completely mine.

clearvars = true;
suffix = '';

if (nargin > 1),
    isvar = true(size(varargin));
    
    i = 1;
    while (i <= length(varargin)),
        if (ischar(varargin{i}) && (varargin{i}(1) == '-')),
            switch lower(varargin{i}),
                case '-noclear',
                    clearvars = false;
                    isvar(i) = false;
                    i = i+1;
                    
                case '-suffix',
                    suffix = varargin{i+1};
                    isvar(i:i+1) = false;
                    i = i+2;
                    
                case '-1',
                    suffix = '1';
                    isvar(i) = false;
                    i = i+1;
                    
                otherwise,
                    error('Unrecognized option %s.',varargin{i});
            end;
        else
            i = i+1;
        end;
    end;
    varargin = varargin(isvar);
end;

if (isempty(varargin)),
    vars = fieldnames(S);
elseif ((length(varargin) == 1) && iscellstr(varargin{1})),
    vars = varargin{1};
else
    vars = varargin;
end;

found = true(size(vars));

if (clearvars),
    vlist = '';
    for i = 1:length(vars),
        vlist = [vlist vars{i} suffix ' '];
    end;
    cmd = ['clear ' vlist(1:end-1) ';'];
    evalin('caller',cmd);
end;

for i = 1:length(vars),
    if (isfield(S,vars{i})),
        assignin('caller',[vars{i} suffix],S.(vars{i}));
    else
        found(i) = false;
    end;
end;

if (nargout == 1),
    varargout{1} = all(found);
elseif (nargout == 2),
    varargout = {all(found) vars(~found)};
elseif ((nargout == 0) && ~all(found)),
    warning('use:missingVars','use did not find all variables requested.');
end;

