function install_here(varargin)
%INSTALL_HERE   Sets up the path structure for Eric Tytell's toolboxes
% install_here(options...)
% Sets up the path structure appropriately for Eric Tytell's toolboxes in
% the current directory.
% Automatically generates a file called "uninstall_here" that removes the
% added paths
%
% Options: 'maxdepth' - maximum depth of subdirectories to include
%          'exclude' - Names of things to exclude.  Uses regexp syntax.
%          Should always have '^\..*' to exclude '.' and '..' on
%          MacOS/Linux systems.

% Mercurial revision hash: $Revision$ $Date$
% See https://bitbucket.org/tytell/matlab/wiki/Home
% Copyright (c) 2011, Eric Tytell <tytell at jhu dot edu>

opt.maxdepth = 2;
opt.exclude = {'^\..*','private','^@'};

%parse varargin without using parsevarargin, since it's not on the path yet
i = 1;
while (i <= length(varargin))
    if (ischar(varargin{i}))
        switch lower(varargin{i})
            case 'maxdepth',
                opt.maxdepth = varargin{i+1};
                i = i+2;
                
            case 'exclude',
                opt.exclude = varargin{i+1};
                if (ischar(opt.exclude))
                    opt.exclude = {opt.exclude};
                end;
                i = i+2;
                
            otherwise
                error('Unrecognized option %s', varargin{i});
        end;
    else
        error('Unrecognized parameter #%d', i);
    end;
end;

basepath = pwd;

butt = questdlg(sprintf('Install at base path %s?',basepath),'Install location','Yes','No','No');
if (~strcmp(butt,'Yes'))
    fprintf('Canceling...\n');
    return;
end;

pathnames = {basepath};

%scan the current directory tree
done = false;
depth = 0;
prevpath = {basepath};
while ((depth <= opt.maxdepth) && ~done)
    done = true;
    nextpath = {};
    for i = 1:length(prevpath),
        files = dir(prevpath{i});
        if (~isempty(files)),
            isdir = cat(1,files.isdir);
            dirs1 = {files(isdir).name};
            dirs1 = cellfun(@(x) (fullfile(prevpath{i},x)), dirs1, 'UniformOutput',false);
        else
            dirs1 = {};
        end;
        
        if (~isempty(dirs1))
            good = true(size(dirs1));
                        
            for j = 1:length(dirs1)
                k = find(dirs1{j} == filesep);
                if (isempty(k))
                    k = 1;
                else
                    k = k(end)+1;
                end;
                for ex = 1:length(opt.exclude)
                    if (regexp(dirs1{j}(k:end),opt.exclude{ex}))
                        good(j) = false;
                    end;
                end;
            end;
            dirs1 = dirs1(good);
                
            if (~isempty(dirs1))
                done = false;
            end;
        end;

        nextpath = [nextpath dirs1];
    end;
    
    pathnames = [pathnames nextpath];
    
    prevpath = nextpath;
    depth = depth + 1;
end;

curpath = regexp(path,pathsep,'split');
newpaths = ~ismember(pathnames,curpath);

if (sum(newpaths) == 0)
    fprintf('No new paths need to be added.\n');
else
    pathnames = pathnames(newpaths);
    
    if (length(pathnames) < 10)
        prompt = [{'Add'},pathnames,{'to the current path?','','(Run uninstall_here to undo)'}];
    else
        extra = sprintf('and %d more',length(pathnames)-9);
        prompt = ['Add',pathnames(1:9),{extra,'to the current path?','','(Run uninstall_here to undo)'}];
    end;
    butt = questdlg(prompt,'Setup path','OK','Cancel','Cancel');
    
    if (strcmp(butt,'OK'))
        try
            fid = fopen('uninstall_here.m','w');
            fprintf(fid,'%% Generated on %s\n', datestr(now));
            for i = 1:length(pathnames)
                fprintf(fid,'fprintf(''Removing ''''%s'''' from the path\\n'');\n', pathnames{i});
                fprintf(fid,'rmpath(''%s'');\n', pathnames{i});
            end;
            fprintf(fid,'savepath;\n');
            fclose(fid);
        catch err
            butt = questdlg({'Could not write to uninstall_here.',...
                'Changes will have to be manually removed.',...
                'Continue?'},...
                'Write error','OK','Cancel','Cancel');
            if ~strcmp(butt, 'OK')
                fprintf('No changes.\n');
                return;
            end
        end
            
        addpath(pathnames{:},'-begin');
        upath = userpath;
        if (~isempty(upath) && ((upath(end) == ':') || (upath(end) == ';')))
            upath = upath(1:end-1);
        end;
        try 
            savepath(fullfile(upath,'pathdef.m'));
        catch err
            fprintf('Error writing path. Writing to startup.m instead.\n');
            fn = fullfile(upath, 'startup.m');
            fid = fopen(fn, 'a');
            fprintf(fid, 'addpath(...\n');
            fprintf(fid, '\t''%s'', ...\n', pathnames{1:end-1});
            fprintf(fid, '\t''%s'', ''-begin'');\n', pathnames{end});
        end
    end;
end;
