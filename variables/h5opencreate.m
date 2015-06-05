function h5opencreate(filename,datasetname,sz, varargin)

opt.Datatype = 'double';
[opt,rest] = parsevarargin(opt,varargin,4,'leaveunknown','typecheck',false);

h5datatypes = struct('H5T_IEEE_F64LE','double', ...
    'H5T_IEEE_F32LE','single', ...
    'H5T_STD_U8LE','uint8', ...
    'H5T_STD_U16LE','uint16', ...
    'H5T_STD_U32LE','uint32', ...
    'H5T_STD_U64LE','uint64', ...
    'H5T_STD_I8LE','int8', ...
    'H5T_STD_I16LE','int16', ...
    'H5T_STD_I32LE','int32', ...
    'H5T_STD_I64LE','int64');

docreate = true;
if exist(filename,'file')
    try
        info = h5info(filename,datasetname);
        if isempty(info.Datatype)
            %it's a group not a dataset
            error('h5opencreate:mismatch','%s is a group, not a dataset',datasetname);
        end
        
        %check for matching datatype
        filedatatype = info.Datatype.Type;
        if isfield(h5datatypes,filedatatype)
            filedatatype = h5datatypes.(filedatatype);
        end
        if ~strcmp(opt.Datatype,filedatatype)
            error('h5opencreate:mismatch','datatypes do not match');
        end
        
        %check for matching size
        if any(~isfinite(sz))
            if ((numel(sz) ~= numel(info.Dataspace.MaxSize)) || ...
                    any(sz ~= info.Dataspace.MaxSize))
                error('h5opencreate:mismatch','size does not match');
            end
        else
            if ((numel(sz) ~= numel(info.Dataspace.Size)) || ...
                    any(sz ~= info.Dataspace.Size))
                error('h5opencreate:mismatch','size does not match');
            end
        end
        docreate = false;
    catch err
        if ~isempty(strfind(err.message, 'name doesn''t exist'))
            % dataset doesn't exist
            docreate = true;
        else
            rethrow(err)
        end
    end
end

if docreate
    h5create(filename,datasetname,sz, 'Datatype',opt.Datatype, ...
        rest{:});
end

        
       
    