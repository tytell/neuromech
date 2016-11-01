function save_data_and_script(scriptname,filename,varargin)

opt.exclude = {};
opt.include = {};
opt = parsevarargin(opt, varargin, 3);

if nargin == 0
    scriptname = '';
    filename = '';
end

if isempty(scriptname) || ~exist(scriptname,'file')
    [fn,pn] = uigetfile('*.m','Script file?');
    if isempty(fn)
        fprint('Canceled. Data not saved.\n');
        return;
    end
    scriptname = fullfile(pn,fn);
end
[pn,fn,ext] = fileparts(scriptname);
if isempty(ext)
    ext = '.m';
end
scriptname = fullfile(pn,[fn ext]);

if isempty(filename)
    [fn,pn] = uiputfile('*.mat','Output file?');
    if isempty(fn)
        fprint('Canceled. Data not saved.\n');
        return;
    end
    filename = fullfile(pn,fn);
end

[pn,fn,ext] = fileparts(filename);
if isempty(pn)
    pn = pwd;
end
filename = fullfile(pn,[fn ext]);

if exist(filename,'file')
    [pn,fn] = fileparts(filename);
    
    [tok,rest] = regexp(fn,'\d+','match','split');
    
    basename = strjoin(rest(1:end-1),tok(1:end-1));
    if length(rest) == length(tok)
        suf = '';
    elseif length(rest) == length(tok)+1
        suf = rest{end};
    else
        error('File name is weird');
    end
    filenames = getfilenames(fullfile(pn,[basename '*' suf '.mat']));
    
    numstr = regexp(filenames{end},[basename '(\d+)' suf '.mat'],'once','tokens');
    numstr = numstr{1};
    num = str2double(numstr) + 1;
    
    filename = sprintf('%s%0*d%s.mat', basename, length(numstr), num, suf);
end

[pn,fn] = fileparts(filename);
scriptsavename = fullfile(pn,[fn '.m']);

if ~isempty(opt.include)
    S = getvar(opt.include{:},'-tostruct');
elseif ~isempty(opt.exclude)
    S = getvar('-all','-except',opt.exclude{:}, '-tostruct');
else
    S = getvar('-all','-tostruct');
end
save(filename,'-struct','S','-v7.3');

copyfile(scriptname,scriptsavename);

fprintf('Data saved to %s and %s\n', filename, scriptsavename);





