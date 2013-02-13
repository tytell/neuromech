function varargout = hg(varargin)
% function hg(parameters...)
% Calls the system hg command with the specified parameters.
%
% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2011, Eric Tytell <tytell at jhu dot edu>


opt.hgcmd = 'hg';
opt.echo = true;
opt.check = false;

[opt,args] = parsevarargin(opt,varargin,'leaveunknown');

if (opt.check)
    [~,result] = system(['which ' opt.hgcmd]);
    if ~isempty(result)
        status = 1;
    else
        status = 0;
    end
else
    if (opt.echo)
        echo = {'-echo'};
    else
        echo = {};
    end;
    
    cmd = sprintf('%s ',opt.hgcmd,args{:},'--config ui.editor=fail');
    
    [status,result] = system(cmd,echo{:});
end

if (nargout == 1)
    varargout = {status};
elseif (nargout == 2)
    varargout = {status,result};
end;

