function S = h5readstruct(filename,varargin)

opt.rootgroup = '/';
opt.structsize = [];
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

if ~exist(filename,'file')
    error('h5readstruct:nofile','File does not exist');
end

info = h5info(filename,opt.rootgroup);

if ~isempty(opt.structsize)
    szhi = opt.structsize;
elseif (~isempty(info.Attributes) && ismember('structsize',{info.Attributes.Name}))
    szhi = h5readatt(filename,opt.rootgroup,'structsize');
else
    szhi = 1;
end
names = {info.Datasets.Name};
szhi = szhi(:)';

C = cell([length(names) szhi]);
for i = 1:length(names)
    nm1 = [opt.rootgroup names{i}];
    D = h5read(filename,nm1);
    
    if (numel(szhi) == 1) && (szhi == 1)
        C{i} = D;
    else
        nd = ndims(D);
        C1 = num2cell(D,1:nd-length(szhi));
        C1 = reshape(C1,[1 szhi]);
        C(i,:) = C1(1,:);
    end
end
S = cell2struct(C,names,1);


