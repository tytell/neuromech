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
        %check to make sure they're all the same size
        nd = cellfun(@ndims,C(i,:));
        nd = max(nd);
        szall = zeros([nd szhi]);
        for j = 1:nd
            szall(j,:) = cellfun(@(x) size(x,j), C(i,:));
        end
        
        sz1 = max(flatten(szall,2:length(szhi)+1),[],2);
        for j = 1:prod(szhi)
            if any(szall(:,j) ~= sz1)
                C1 = C{i,j};
                for k = 1:nd
                    C1(end+1:sz1(k),:) = NaN;
                    C1 = shiftdim(C1,1);
                end
                C{i,j} = C1;
            end
        end

        sz1 = sz1';
        D = cat(length(sz1)+1,C{i,:});
        D = reshape(D,[sz1 szhi]);
        
        if ~isreal(D)
            D = shiftdim(D,-1);
            D(2,:) = imag(D(1,:));
            D(1,:) = real(D(1,:));
            sz1 = [2 sz1];
            Disreal = 0;
        else
            Disreal = 1;
        end
        nm1 = [opt.rootgroup names{i}];
        if ~ismember(names{i},existnames)
            h5create(filename,nm1, [sz1 Inf(size(szhi))], 'ChunkSize',...
                [sz1 ones(size(szhi))], 'Datatype',class(D));
        end
        h5write(filename,nm1, D, ones(1,length(sz1)+length(szhi)), [sz1 szhi]);
        h5writeatt(filename,nm1,'isreal',Disreal);
    else
        warning('h5writestruct:type','Skipping field %s because type is not numeric',names{i});
    end
end
h5writeatt(filename,opt.rootgroup,'structsize',szhi);

        
        