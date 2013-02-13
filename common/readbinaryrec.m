function S = readbinaryrec(fid, rec, N, varargin)

opt.order = 'record';

opt = parsevarargin(opt,varargin,4);

C = struct2cell(rec);
varnames = fieldnames(rec);

nvar = length(C);
varlen = zeros(length(C),1);
varnum = zeros(length(C),1);
vartype = cell(length(C),1);
variablelenref = cell(length(C),1);

maxdims = max(cellfun(@ndims,C));
matsz = ones(length(C),maxdims);

for i = 1:length(C),
    [varlen(i),varnum(i),vartype{i}, ref1] = bytelen(C{i});
    matsz(i,1:ndims(C{i})) = size(C{i});
    
    if (varnum(i) == -1),
        refind = -1 * ones(1,length(ref1));
        for j = 1:length(ref1),
            if (~isempty(ref1{j}))
                refind1 = find(strcmp(ref1{j}, varnames));
                
                if (isempty(refind1))
                    error('readbinaryrec:badlengthspec',...
                        'Length specifier field %s does not exist',ref1{j});
                end;
                refind(j) = refind1;
            end;
        end;
        variablelenref{i} = refind;
    end;
end;

S = cell(length(C),N);
for j = 1:N,
    for i = 1:nvar,
        if (varnum(i) > 0)
            if (strcmp(vartype{i},'char'))
                S{i,j} = char(zeros(matsz(i,:)));
            else
                S{i,j} = zeros(matsz(i,:),vartype{i});
            end;
        end;
    end;
end;
    
if (all(varnum > 0))
    %easier to read fixed length records...
    reclen = sum(varlen.*varnum);
    
    switch lower(opt.order),
        case 'record',
            for j = 1:N,
                buf = fread(fid, reclen, '*uint8');

                for i = 1:nvar,
                    S{i,j} = typecast(buf(k+(1:varlen(i)*varnum(i))),vartype{i});
                    k = k+varlen(i)*varnum(i);
                end;
                        
            end;
            
        case 'value',
            startpos = ftell(fid);
            for i = 1:nvar,
                spec = sprintf('%d*%s',varnum(i),vartype{i});
                S1 = fread(fid, N, spec, reclen - varlen(i)*varnum(i));
                S(i,:) = mat2cell(S1,1,varlen(i)*ones(1,N));
                
                startpos = startpos + varlen(i)*varnum(i);
                fseek(fid, startpos, 'bof');
            end;
    end;
else
    if (strcmp(opt.order,'value'))
        error('readbinaryrec:valuepriority','Cannot read in order of value with variable length records');
    end;
    
    for i = 1:nvar,
        for j = 1:N,
            if (varnum(i) > 0)
                S{i,j}(:) = fread(fid, varnum(i), ['*' vartype{i}]);
            else
                ind = variablelenref{i};
                lens = -1*ones(size(ind));
                good = ind > 0;
                speclen = cat(1,S{ind(good),j});
                lens(good) = speclen;
                
                S{i,j} = readvariablelen(fid, varlen(i),vartype{i}, lens);
            end;
        end;
    end;
end;

S = cell2struct(S,varnames,1);



function [len,N,tp,ref] = bytelen(val)

ref = {};
tp = class(val);
switch tp,
    case {'int8','uint8'},
        len = 1;
        N = numel(val);
    case {'int16','uint16'},
        len = 2;
        N = numel(val);
    case {'int32','uint32'},
        len = 4;
        N = numel(val);
    case {'int64','uint64'},
        len = 8;
        N = numel(val);
    case 'double',
        len = 8;
        N = numel(val);
    case 'single',
        len = 4;
        N = numel(val);
    case 'char',
        len = 1;
        N = numel(val);
        if (N == 0)
            %null terminated character string
            N = -1;
            ref = '';
        end;
        
    case 'cell',
        [len, N1, tp, ref1] = bytelen(val{1});
        N = -1;
        
        switch lower(val{2}),
            case 'nullterminated',
                ref2 = '';
            case 'specified',
                ref2 = val{3};
        end;
        
        if (N1 == -1)
            %handle nested variable length things...
            ref = [ref2 ref1];
        else
            ref = {ref2};
        end;
        
    otherwise,
        error('Unsupported data type %s',tp);
end;



function val = readvariablelen(fid, itemlen,itemtype, lens)

if ((length(lens) > 1) && all(lens > 0)),
    %all lengths specified
    val = fread(fid,prod(lens),['*' itemtype]);
    val = reshape(val,lens);
elseif (length(lens) == 1)
    if (lens == -1)
        %null terminated
        len = 0;
        if (strcmp(itemtype,'char'))
            zero = char(0);
        else
            zero = zeros(1,1,itemtype);
        end;
        while (fread(fid,1,['*' itemtype]) ~= zero)
            len = len + 1;
        end;
        fseek(fid,-(len+1)*itemlen,'cof');
        val = fread(fid, len, ['*' itemtype]);
        fseek(fid,itemlen,'cof');     % skip over null terminator
    else
        val = fread(fid, lens, ['*' itemtype]);
    end;
else
    if (lens(1) == -1)
        if (strcmp(itemtype,'char'))
            zero = char(0);
        else
            zero = zeros(1,1,itemtype);
        end;

        val = cell(1,0);
        done = false;
        a = 1;
        while (~done)
            val{a} = readvariablelen(fid, itemlen,itemtype, lens(2:end));
            isstop = fread(fid,1,['*' itemtype]);
            if (isstop == zero)
                done = true;
            else
                fseek(fid,-itemlen,'cof');
            end;
        end;
    else
        val = cell(1,lens(1));
        for i = 1:lens(1),
            val{i} = readvariablelen(fid, itemlen,itemtype, lens(2:end));
        end;
    end;
end;


