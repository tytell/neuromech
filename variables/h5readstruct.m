function S = h5readstruct(filename,varargin)

opt.rootgroup = '/';
opt = parsevarargin(opt,varargin,3);

if (isempty(opt.rootgroup))
    opt.rootgroup = '/';
else
    if (opt.rootgroup(1) ~= '/')
        opt.rootgroup = ['/' opt.rootgroup];
    end
    if (opt.rootgroup(end) ~= '/')
        opt.rootgroup = [opt.rootgroup '/'];
    end
end

info = h5info(filename,opt.rootgroup);
names = {info.Datasets.Name};

C = cell(length(names),1);
for i = 1:length(names)
    nm1 = [opt.rootgroup names{i}];
    C{i} = h5read(filename,nm1);
end
S = cell2struct(C,names,1);


