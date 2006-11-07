function checkfcnshadowing(varargin)
% function checkfcnshadowing(files/dir)
% Checks for shadowed functions.  Pass in a directory name or a list of
% file names.  With no arguments checks the current directory.  Gives a
% list of functions that shadow others, assuming that you can figure out
% which ones are shadowed.
% Option --verbose displays exactly which functions are shadowed.

ind = strmatch('--verbose',varargin,'exact');
if (~isempty(ind)),
    isverbose = true;
    varargin = varargin([1:ind-1 ind+1:end]);
else
    isverbose = false;
end;

files = {};
if (nargin == 0),
    dir = '.';
    fprintf('Checking current directory (%s):\n',pwd);
elseif ((nargin == 1) && exist(varargin{1},'dir')),
    dir = varargin{1};
    fprintf('Directory %s:',dir);
else
    files = varargin{1};
end;

if (isempty(files)),
    files = getfilenames(fullfile(dir,'*.m'));
end;

for f = 1:length(files),
    [path,name] = fileparts(files{f});
    shadow = which('-all',name);

    if (isverbose),
        if (length(shadow) > 1),
            fprintf('%s - Shadows\n', name);
            fprintf('     %s\n',shadow{2:end});
        else
            fprintf('%s - OK\n',files{f});
        end;
    elseif (length(shadow) > 1),
        fprintf('%s - Shadows %d functions.\n',name,length(shadow)-1);
    end;
end;
fprintf('Done.\n');


            
