function hg(varargin)
% function hg(parameters...)
% Calls the system hg command with the specified parameters.
%
% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2010, Eric Tytell <tytell at jhu dot edu>


opt.hgpath = '/usr/local/bin/hg';

[opt,args] = parsevarargin(opt,varargin,'leaveunknown');

cmd = sprintf('%s ',opt.hgpath,args{:});

system(cmd);
