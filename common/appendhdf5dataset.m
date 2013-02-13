function id = appendhdf5dataset(filename,datasetname, data,dim, varargin)

opt.extendable = true;
opt.chunksize = [];
opt.truncate = false;
opt.ndims = [];
opt.fillvalue = NaN;

if (isstruct(filename))
    id = filename;
    init = false;
else
    id = struct('file',[]);
    init = true;
end;
if (isempty(datasetname))
    idsetname = 'set';
else
    idsetname = datasetname;
    idsetname(idsetname == '/') = '_';
    if (idsetname(1) == '_')
        idsetname = idsetname(2:end);
    end;
end;
if ((nargin == 3) || ((nargin > 3) && isempty(dim)))
    dim = ndims(data)+1;
end;

if (ischar(datasetname) && ismember(lower(datasetname),{'close'}))
    cmd = lower(datasetname);
    switch cmd,
        case 'close',
            fn = fieldnames(id);
            for i = 1:length(fn),
                if (~strcmp(fn{i},'file'))
                    H5D.close(id.(fn{i}));
                end;
            end;
            H5F.close(id.file);
            id = struct([]);
            return;
    end;
end;

opt = parsevarargin(opt,varargin,5);

if (~isfield(id,idsetname))
    init = true;
end;

if (init)
    if (isempty(id.file))
        if (exist(filename,'file'))
            id.file = H5F.open(filename, 'H5F_ACC_RDWR', 'H5P_DEFAULT');
        else
            if (opt.truncate)
                createmode = 'H5F_ACC_TRUNC';
            else
                createmode = 'H5F_ACC_EXCL';
            end;
            id.file = H5F.create(filename, createmode,'H5P_DEFAULT','H5P_DEFAULT');
        end;
    end;
    
    %try to open the dataset
    [id.(idsetname), err] = hdf5err(@H5D.open, id.file, datasetname);
    if (~isempty(err))
        %if we couldn't, then create it
        
        %set the base size to the size of data
        sz = size(data);
        
        %make sure we have enough dimensions
        if (~isempty(opt.ndims))
            nd = opt.ndims;
        elseif (dim > ndims(data))
            nd = dim;
        else
            nd = ndims(data);
        end;
        sz(end+1:nd) = 1;
        
        %start with 0 size along the dimension we're appending
        sz(dim) = 0;
        
        %set the maxsz appropriately
        maxsz = sz;
        if (length(opt.extendable) == nd)
            maxsz(opt.extendable) = H5ML.get_constant_value('H5S_UNLIMITED');
        elseif ((numel(opt.extendable) == 1) && opt.extendable)
            maxsz(:) = H5ML.get_constant_value('H5S_UNLIMITED');
        end;
        maxsz(dim) = H5ML.get_constant_value('H5S_UNLIMITED');
        
        %and create the dataset
        filespace = H5S.create_simple(nd, sz, maxsz);
        
        fillval = opt.fillvalue;
        switch class(data),
            case 'uint8',
                htype = H5T.copy('H5T_NATIVE_B8');
                fillval = uint8(0);
            case 'uint16',
                htype = H5T.copy('H5T_NATIVE_B16');
                fillval = uint16(0);
            case 'uint32',
                htype = H5T.copy('H5T_NATIVE_B32');
                fillval = uint32(0);
            case 'double',
                htype = H5T.copy('H5T_NATIVE_DOUBLE');
            otherwise
                error('Unrecognized type for the data');
        end;
        
        if (any(opt.extendable) && isempty(opt.chunksize))
            chunksize = size(data);
            chunksize(end+1:nd) = 1;
        else
            chunksize = opt.chunksize;
        end;
        
        cparms = H5P.create('H5P_DATASET_CREATE');
        H5P.set_chunk(cparms,chunksize);
        H5P.set_fill_value(cparms,htype,fillval);
        
        id.(idsetname) = H5D.create(id.file, datasetname, htype,filespace, cparms);
        H5T.close(htype);
    end;
end;

%check the size of the dataset
filespace = H5D.get_space(id.(idsetname));
[nd,datasetsz] = H5S.get_simple_extent_dims(filespace);

sz = size(data);
if (ndims(data) < nd)
    sz(end+1:nd) = 1;
end;

%extend along dim
start = zeros(1,nd);
start(dim) = datasetsz(dim);
datasetsz(dim) = datasetsz(dim) + sz(dim);

%and check for any other dimensions that are too small
small = datasetsz < sz;
if (any(small))
    datasetsz(small) = sz(small);
end;
    
%do the extension
H5D.set_extent(id.(idsetname), datasetsz);
    
%reopen the dataspace
H5S.close(filespace);
filespace = H5D.get_space(id.(idsetname));

%select the hyperslab and write
memspace = H5S.create_simple(nd, sz,sz);
H5S.select_hyperslab(filespace, 'H5S_SELECT_SET',start,[],sz,[]);

data = permute(data,[2 1 3:nd]);

H5D.write(id.(idsetname),'H5ML_DEFAULT',memspace, filespace, 'H5P_DEFAULT', data);

%close necessary things
H5S.close(filespace);
H5S.close(memspace);

if (nargout == 0)
    H5D.close(id.(idsetname));
    H5F.close(id.file);
end;
