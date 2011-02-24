function packagefcn(fcns,dest, varargin)
% function packagefcn(fcns,dest, options...)
% Finds the dependencies of a particular function(s) and copies them all
% into a separate directory.
% Option 'maintaindirstruct' keeps the overall directory structure, but
% only copies the relevant files.

% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2010, Eric Tytell <tytell at jhu dot edu>

opt.maintaindirstruct = true;
opt.dryrun = false;

opt = parsevarargin(opt, varargin, 2);

if ((nargin == 0) || isempty(fcns)),
    [fcns,~] = uigetfile('*.m','Select function(s) to package', 'MultiSelect','on');
    if (~ischar(fcns) && (fcns == 0))
        return;
    end;
    
    dest = uigetdir('','Select destination direction');
    if (~ischar(dest) && (dest == 0))
        return;
    end;
end;

if (~iscell(fcns)),
    fcns = {fcns};
end;

package = {};
for i = 1:length(fcns),
    fullname = which(fcns{i});
    [pn,fn,~] = fileparts(fullname);
    cd(pn);
    
    dep = depfun(fn);
    
    builtin = strncmp(matlabroot,dep,length(matlabroot));
    isbuiltin = false(size(dep));
    isbuiltin(builtin) = true;
    
    package = union(package,dep(~isbuiltin));
end;

if (opt.maintaindirstruct),
    makedirs = {};
    for i = 1:length(package),
        tok = regexp(package{i},'/', 'split');
        tok = tok(2:end-1);

        j = 1;
        ismatch = false;
        while ((j <= length(makedirs)) && ~ismatch),
            if (length(tok) == length(makedirs{j})),
                ismatch = all(cellfun(@strcmpi, tok,makedirs{j}));
            else
                ismatch = false;
            end;
            j = j+1;
        end;
        
        if (~ismatch),
            makedirs{end+1} = tok;      %#ok
        end;
    end;
    
    len = cellfun(@length,makedirs);
    [len,ord] = sort(len);
    makedirs = makedirs(ord);
    
    minlen = len(1);
    common = true(1,minlen);
    i = 1;
    while ((i <= minlen) && common(i)),
        j = 2;
        while ((j <= length(makedirs)) && common(i)),
            common(i) = common(i) & strcmpi(makedirs{j}{i},makedirs{1}{i});
            j = j+1;
        end;
        i = i+1;
    end;
    
    common = last(common);
    commonbasedir = fullfile('/',makedirs{1}{1:common});
    
    for i = 1:length(makedirs),
        if (length(makedirs{i}) > common),
            dirname = fullfile(dest,makedirs{i}{common+1:end});
        
            if (opt.dryrun)
                fprintf('mkdir(%s)\n',dirname);
            else
                mkdir(dirname);
            end;
        end;
    end;
    
    cutlen = length(commonbasedir);
    for i = 1:length(package),
        [pn,~,~] = fileparts(package{i});
        pn = pn(cutlen+1:end);
        
        if (opt.dryrun)
            fprintf('copyfile(%s,%s)\n', package{i}, fullfile(dest,pn));
        else
            copyfile(package{i},fullfile(dest,pn));
        end;
    end;
else                       
    for i = 1:length(package),
        if (opt.dryrun)
            fprintf('copyfile(%s,%s)\n', package{i}, dest);
        else
            copyfile(package{i},dest);
        end;
    end;
end;



