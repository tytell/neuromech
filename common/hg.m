function varargout = hg(varargin)
% function hg(parameters...)
% Calls the system hg command with the specified parameters.
%
% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2011, Eric Tytell <tytell at jhu dot edu>


opt.hgpath = '/usr/local/bin/hg';
opt.echo = true;

[opt,args] = parsevarargin(opt,varargin,'leaveunknown');

if (opt.echo)
    echo = {'-echo'};
else
    echo = {};
end;
    
cmd = sprintf('%s ',opt.hgpath,args{:},'--config ui.editor=fail');

[status,result] = system(cmd,echo{:});

if (nargout == 1)
    varargout = {status};
elseif (nargout == 2)
    varargout = {status,result};
end;

