function hg(varargin)

opt.hgpath = '/usr/bin/hg';

[opt,args] = parsevarargin(opt,varargin,'leaveunknown');

cmd = sprintf('%s ',opt.hgpath,args{:});

system(cmd);
