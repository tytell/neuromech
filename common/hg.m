function hg(varargin)

opt.hgpath = '/usr/local/bin/hg';

[opt,args] = parsevarargin(opt,varargin,'leaveunknown');

cmd = sprintf('%s ',opt.hgpath,args{:});

system(cmd);
