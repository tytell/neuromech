function packagefcn(fcns,dest, varargin)
% function packagefcn(fcns,dest, options...)

opt.maintaindirstruct = false;

i = 1;
while (i <= length(varargin)),
    switch lower(varargin{i}),
      case {'maintaindirstruct'},
        opt.(varargin{i}) = true;
        i = i+1;
        
      otherwise,
        error('Unrecognized option %s', varargin{i});
    end;
end;

if (nargin == 0),
    [fcns,pn] = uigetfile('*.m','Select function(s) to package', 'MultiSelect','on');
    
    dest = uigetdir('','Select destination direction');
end;

if (~iscell(fcns)),
    fcns = {fcns};
end;

package = {};
for i = 1:length(fcns),
    fullname = which(fcns{i});
    [pn,fn,ext] = fileparts(fullname);
    cd(pn);
    
    dep = depfun(fn);
    
    builtin = strmatch(matlabroot,dep);
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
            makedirs{end+1} = tok;
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
        
            mkdir(dirname);
        end;
    end;
    
    cutlen = length(commonbasedir);
    for i = 1:length(package),
        [pn,fn,ext] = fileparts(package{i});
        pn = pn(cutlen+1:end);
        
        copyfile(package{i},fullfile(dest,pn));
    end;
else                       
    for i = 1:length(package),
        copyfile(package{i},dest);
    end;
end;



