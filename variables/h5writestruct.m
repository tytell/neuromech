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

for i = 1:size(C,1)
    if isnumeric(C{i,1})
        sz1 = size(C{i,1});
        D = cat(length(sz1)+1,C{i,:});
        D = reshape(D,[sz1 szhi]);
        
        nm1 = [opt.rootgroup names{i}];
        h5create(filename,nm1, size(D), 'Datatype',class(D));
        h5write(filename,nm1, D);
    else
        warning('h5writestruct:type','Skipping field %s because type is not numeric',names{i});
    end
end
        
        