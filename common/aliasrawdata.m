function aliasrawdata
% function aliasrawdata
% Creates symbolic links to data files on an external hard drive
%
% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2010, Eric Tytell <tytell at jhu dot edu>

if (~ismac && ~isunix),
    error('Aliases do not work well under Windows');
end;

[fn,pn] = uigetfile('*','Select raw data files','/Volumes', ...
    'MultiSelect','on');
if (~iscell(fn)),
    if (ischar(fn)),
        fn = {fn};
    elseif (isnumeric(fn) && (fn == 0)),
        return;
    end;
end;

outdir = uigetdir('','Select output directory');
if (~ischar(outdir) && (outdir == 0)),               % canceled
    return;
end;

if (outdir(end) ~= '/'),
    outdir(end+1) = '/';
end;

overwriteall = false;
for i = 1:length(fn),
    fn1 = fullfile(pn,fn{i});
    
    cmd = ['ln -sf ' fn1 ' ' outdir];
    if (~overwriteall && exist(fullfile(outdir,fn{i}),'file')),
        butt = questdlg(sprintf('Alias to %s exists.  Overwrite?',fn{i}), ...
            'Overwrite?','This file only','All files','Skip','Skip');
        
        switch butt,
            case 'All files',
                overwriteall = true;
            case 'Skip',
                cmd = '';
        end;
    end;
        
    if (~isempty(cmd)),
        system(cmd);
    end;
end;

    