function dirs = getdirectorynames(d)

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
