function dirs = getdirectorynames(d)
% function dirs = getdirectorynames(d)
% Returns directories on a given path (potentially with wildcards).
%
% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2010, Eric Tytell <tytell at jhu dot edu>

[pathnm,~] = fileparts(d);

files = dir(d);
if (~isempty(files)),
    isdir = cat(1,files.isdir);
    dirs = {files(isdir).name};
    dirs = cellfun(@(x) (fullfile(pathnm,x)), dirs, 'UniformOutput',false);
else
	dirs = {};
end;

dirs = sortfiles(dirs);
