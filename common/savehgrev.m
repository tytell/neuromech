function out = savehgrev(out, varargin)
%SAVEHGREV  Produces an HGREV structure with the current changeset
%   savehgrev(outfile, options...)
%     or
%   outstruct = savehgrev(outstruct, options...)
%
%  The HGREV structure contains the revision number, the changeset hash,
%  and the description, branch, and date information.  Useful to replicate
%  calculations that produced a particular data file.
%
% OPTIONS
% General options:
%   'caller' - Identifies the calling function that generates a data file.
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
%   HG

% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2011, Eric Tytell <tytell at jhu dot edu>

opt.uncommitted = 'warning';
opt.untracked = 'none';
opt.nottip = 'warning';
opt.caller = '';

opt = parsevarargin(opt,varargin, 2);

[s,res] = hg('summary','echo',false);

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
    if (~isempty(opt.caller))
        [~,stat] = hg('status -A',opt.caller,'echo',false);
        sttok = regexp(stat, '^([MARC!?I])', 'once','tokens');
   
        if (sttok{1} ~= 'C')
            isuncom = true;
        end;
    else
        [~,stat] = hg('status -mard','echo',false);
        if (~isempty(stat))
            isuncom = true;
        end;
    end;
    if (isuncom)
        feval(warnfunc,'savehgrev:uncommitted','Files are uncommitted in the current directory.');
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
    [~,stat] = hg('status -u','echo',false);
    if (~isempty(stat))
        feval(warnfunc,'savehgrev:untracked','Files are untracked in the current directory.');
    end;
end;
   
%get changeset info
cstok = regexp(res, 'parent: (\d+):([0-9a-f]+)\s(\w+)?', 'tokens','once');
currev = str2double(cstok{1});
%get the full changeset hash
[~,cshash] = hg('log -r ',cstok{1},' --template ''{node}''','echo',false);

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
        feval(warnfunc,'savehgrev:nottip','Current directory is not at tip revision');
    end;
end;

[~,branch] = hg('log -r ',cstok{1},' --template ''{branches}''','echo',false);
[~,desc] = hg('log -r ',cstok{1},' --template ''{desc}''','echo',false);
[~,revdate] = hg('log -r ',cstok{1},' --template ''{date|isodate}''','echo',false);

HGREV.rev = currev;
HGREV.changeset = cshash;
HGREV.description = desc;
HGREV.branch = branch;
HGREV.date = revdate;
if (~isempty(opt.caller))
    HGREV.caller = opt.caller;
end;

if (ischar(out) && exists(out,'file'))
    save(out,'HGREV','-append');
elseif (isstruct(out))
    out.HGREV = HGREV;
end;








