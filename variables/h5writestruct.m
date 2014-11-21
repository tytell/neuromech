function h5writestruct(filename,S, varargin)

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

szhi = size(S);
C = struct2cell(S);
names = fieldnames(S);

if isempty(S)
    return;
end
    
if exist(filename,'file')
    try
        info = h5info(filename,opt.rootgroup);
        existnames = {info.Datasets.Name};
    catch me
        existnames = {};
    end
else
    existnames = {};
end
        
for i = 1:size(C,1)
    if isnumeric(C{i,1})
        sz1 = size(C{i,1});
        D = cat(length(sz1)+1,C{i,:});
        D = reshape(D,[sz1 szhi]);
        
        nm1 = [opt.rootgroup names{i}];
        if ~ismember(names{i},existnames)
            h5create(filename,nm1, [sz1 Inf(size(szhi))], 'ChunkSize',...
                [sz1 ones(size(szhi))], 'Datatype',class(D));
        end
        h5write(filename,nm1, D, ones(1,length(sz1)+length(szhi)), [sz1 szhi]);
    else
        warning('h5writestruct:type','Skipping field %s because type is not numeric',names{i});
    end
end
h5writeatt(filename,opt.rootgroup,'structsize',szhi);

        
        