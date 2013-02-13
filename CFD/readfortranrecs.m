function S = readfortranrecs(filename, datatype, varargin)
% function S = readfortranrecs(filename, datatype, varargin)
% Reads Fortran binary "rec" output.  The format of the record is specified
% in datatype, which is a struct containing elements of the appropriate
% size.
%
% datatype = struct('A',uint8(0),'B',uint16(0),'C',uint8(0),
% 'D',double(0));
% S = readfortranrecs('test.rec',datatype);
%
% That would read a record consisting of a uint8, a uint16, another uint8,
% and a double.  Note that the order of fields in the structure matters.

opt.method = 'valskip';

opt = parsevarargin(opt,varargin,2);

fid = fopen(filename,'r');
fseek(fid, 0, 'eof');
filesize = ftell(fid);
fseek(fid, 0, 'bof');

C = struct2cell(datatype);
nvar = length(C);
varlen = zeros(length(C),1);
vartype = cell(length(C),1);
for i = 1:length(C),
    vartype{i} = class(C{i});
    switch vartype{i},
        case {'int8','uint8'},
            varlen(i) = 1;
        case {'int16','uint16'},
            varlen(i) = 2;
        case {'int32','uint32'},
            varlen(i) = 4;
        case {'int64','uint64'},
            varlen(i) = 8;
        case 'double',
            varlen(i) = 8;
        case 'single',
            varlen(i) = 4;
        otherwise,
            error('Unsupported data type %s',vartype{i});
    end;
end;
estreclen = sum(varlen);

nrecs = filesize/(estreclen + 8);
if (nrecs ~= floor(nrecs))
    warning('readfortranrecs:trunc',...
        'Records are truncated or variable length.  Attempting to read');
    nrecs = floor(nrecs);
end;

S = cell(length(C),1);
for i = 1:length(C),
    S{i} = zeros(1,nrecs,class(C{i}));
end;
    
switch lower(opt.method),
    case 'record',
        timedWaitBar(0,'Loading records');
        
        j = 1;
        done = false;
        while (~done),
            buf = fread(fid, estreclen+8, '*uint8');
            if (isempty(buf)),
                done = true;
            else
                recl = typecast(buf(1:4),'uint32');
                if (recl == estreclen),
                    k = 4;
                    for i = 1:nvar,
                        S{i}(j) = typecast(buf(k+(1:varlen(i))),class(C{i}));
                        k = k+varlen(i);
                    end;
                    
                    chk = typecast(buf(k+1:end),'uint32');
                    if (chk ~= recl)
                        warning('File corrupted at record %d.  Ending...',j);
                        done = true;
                    end;
                else
                    warning('Bad record length at record %d.  Ending...',j);
                    done = true;
                end;
                j = j+1;
            end;
            if (mod(j,100) == 0)
                if (~timedWaitBar(j/nrecs))
                    done = true;
                end;
            end;
        end;
        timedWaitBar(1);
    case 'valskip',
        totalreclen = estreclen+8;
        
        recl = fread(fid, inf, '*uint32', totalreclen - 4);
        
        if (all(recl == estreclen)),
            startpos = 4;
            for i = 1:nvar,
                fseek(fid, startpos, 'bof');
                
                S{i} = fread(fid, inf, ['*' class(C{i})], totalreclen - varlen(i));
                
                startpos = startpos + varlen(i);
            end;
        end;
end;

fclose(fid);

S = cell2struct(S,fieldnames(datatype),1);

    