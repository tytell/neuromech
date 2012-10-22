function out = savehgrev(out, varargin)
%SAVEHGREV  Produces an HGREV structure with the current changeset
%   savehgrev(outfile, [repo], options...)
%     or
%   outstruct = savehgrev(outstruct, [repo], options...)
%
%  The HGREV structure contains the revision number, the changeset hash,
%  and the description, branch, and date information.  Useful to replicate
%  calculations that produced a particular data file.  Takes an optional
%  parameter, repo, that specifies the repository directory.
%
%  If the datafile specified in the options is a .mat file and has a 
%  structure called HGREV, we save that in HGREV.sub.  So if you have
%  several layers of processing, you may end up with several levels of
%  HGREV.sub.sub....sub.
%
% OPTIONS
% General options:
%   'datafile' - Identifies the input data file to the calling function
%     (ie, the raw data).
%   'savemainrepo' - Saves the state of a repository containing common
%     functions.
%   'mainrepolocation' - Location of the common repository.  Default:
%     '~/matlab'.
%   'disablewarnings' - Disable all the warning options
%   'hash' - Save the SHA hash of the data file.  If false, just save the
%     size and date modified.
%
% Warning options:  Each of these can have the value 'error','warning', or
%      'none', and produces an error, warning, or nothing.
%   
%   'uncommitted' - Warn if the calling function is uncommitted (or if any
%     files are uncommitted, if the caller isn't specified).  Default:
%     'warn'
%   'untracked' - Warn if files are untracked.  Default: 'none'
%   'nottip' - Warn if the current directory is not at the tip.  Default:
%     'warn'
%
% SEE ALSO
%   HG, VERIFYHGREV

% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2011, Eric Tytell <tytell at jhu dot edu>

opt.uncommitted = 'warning';
opt.untracked = 'none';
opt.nottip = 'warning';
opt.datafile = '';
opt.savemainrepo = true;
opt.mainrepolocation = '~/matlab';
opt.disablewarnings = false;
opt.hash = true;
opt.note = '';

if ~hg('check')
    warning('savehgrev:nohg','Mercurial not found');
    out = struct([]);
    return;
end

if ((length(varargin) >= 1) && exist(varargin{1},'dir'))
    repo = {'-R ',varargin{1}};
    repoloc = varargin{1};
    p = 2;
else
    repo = {};
    repoloc = 'current directory';
    p = 1;
end;

opt = parsevarargin(opt,varargin(p:end), p+1, 'typecheck',false);
if (opt.disablewarnings)
    opt.uncommitted = 'none';
    opt.untracked = 'none';
    opt.nottip = 'none';
end;

st = dbstack(1);
callfcns = {st.file};
    
[s,res] = hg(repo{:},'summary','echo',false);
if (s ~= 0)
    error('savehgrev:hgerror','Error running hg');
end;

%check for uncommitted
switch lower(opt.uncommitted)
    case 'warning',
        warnfunc = @warning;
    case 'error',
        warnfunc = @error;
    case 'none',
        warnfunc = [];
    otherwise,
        error('savehgrev:badoption','Unrecognized value ''%s'' for uncommitted option', opt.uncommitted);
end;
if (~isempty(warnfunc))
    isuncom = false;
    [~,stat] = hg(repo{:},'status -mard','echo',false);
    if (~isempty(stat))
        feval(warnfunc,'savehgrev:uncommitted','Files are uncommitted in %s.',repoloc);
    end;
end;

%check for untracked
switch lower(opt.untracked)
    case 'warning',
        warnfunc = @warning;
    case 'error',
        warnfunc = @error;
    case 'none',
        warnfunc = [];
    otherwise,
        error('savehgrev:badoption','Unrecognized value ''%s'' for untracked option', opt.untracked);
end;
if (~isempty(warnfunc))
    [~,stat] = hg(repo{:},'status -u','echo',false);
    if (~isempty(stat))
        feval(warnfunc,'savehgrev:untracked','Files are untracked in %s.',repoloc);
    end;
end;
   
%get changeset info
cstok = regexp(res, 'parent: (\d+):([0-9a-f]+)\s(\w+)?', 'tokens','once');
currev = str2double(cstok{1});
%get the full changeset hash
[~,cshash] = hg(repo{:},'log -r ',cstok{1},' --template ''{node}''','echo',false);

%check to see if we're not at the tip
if (~strcmp(cstok{3},'tip'))
    switch lower(opt.nottip)
        case 'warning',
            warnfunc = @warning;
        case 'error',
            warnfunc = @error;
        case 'none',
            warnfunc = [];
        otherwise,
            error('savehgrev:badoption','Unrecognized value ''%s'' for ''nottip'' option', opt.untracked);
    end;
    if (~isempty(warnfunc))
        feval(warnfunc,'savehgrev:nottip','%s is not at tip revision',repoloc);
    end;
end;

[~,branch] = hg(repo{:},'log -r ',cstok{1},' --template ''{branches}''','echo',false);
[~,desc] = hg(repo{:},'log -r ',cstok{1},' --template ''{desc}''','echo',false);
[~,revdate] = hg(repo{:},'log -r ',cstok{1},' --template ''{date|isodate}''','echo',false);

HGREV.rev = currev;
HGREV.changeset = cshash;
HGREV.description = desc;
HGREV.branch = branch;
HGREV.date = revdate;
HGREV.caller = callfcns;
HGREV.note = opt.note;

if (~isempty(opt.datafile))
    if (ischar(opt.datafile))
        datafiles = {opt.datafile};
    else
        datafiles = opt.datafile;
    end;
    issub = false;
    sub = cell(1,length(datafiles));
    for i = 1:length(datafiles),
        dfinfo1 = dir(datafiles{i});
        datafile1 = datafiles{i};
        if (isempty(dfinfo1))
            datafile1 = [datafile1 '.mat'];
            dfinfo1 = dir(datafile1);
            if (isempty(dfinfo1))
                error('savehgrev:unknowndatafile','Data file %s not found',datafiles{i});
            end;
            ismatfile = true;
        else
            ismatfile = strcmp(datafile1(end-3:end),'.mat');
        end;
        %the output of dir doesn't include the path, so add it on here
        dfinfo1.name = datafiles{i};

        if (opt.hash)
            %get SHA1 hash of data file
            [s,dfhash1] = system(['openssl sha1 ' datafile1]);
            if (s == 0)
                hash1 = regexp(dfhash1,'=\s*([0-9a-f]+)','tokens','once');
                dfinfo1.sha1 = hash1{1};
            end;
        end;
        
        %also check to see if it has an HGREV
        if (ismatfile && ~isempty(who('-file',datafile1,'HGREV')))
            F = load(datafile1,'HGREV');
            sub{i} = F.HGREV;
            issub = true;
        end;
        dfinfo{i} = dfinfo1;
    end;
    dfinfo = makestructarray(dfinfo{:});
    
    if (issub)
        HGREV.sub = makestructarray(sub{:});
    end;
    
    HGREV.datafiles = dfinfo;
end;

if (opt.savemainrepo)
    opt1 = opt;
    opt1.savemainrepo = false;
    opt1.datafile = '';
    kv = opt2keyvalue(opt1);
    
    mainrepo = savehgrev([], opt.mainrepolocation, kv{:});
    mainrepo.caller = {};
    
    HGREV.main = mainrepo;
    HGREV.main.location = opt.mainrepolocation;
end;

if (ischar(out) && exist(out,'file'))
    if (~isempty(who('HGREV','-file',out)))
        F = load(out,'HGREV');
        F.HGREV = makestructarray(F.HGREV,HGREV);
    else
        F.HGREV = HGREV;
    end;
    save(out,'-struct','F','-append');
elseif (isempty(out) || isstruct(out))
    out = HGREV;
end;








