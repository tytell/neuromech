function data = readtecplot(fn, varargin)

opt.binary = false;

opt = parsevarargin(opt,varargin,2);

if (opt.binary)
    tecplotstring = {int32(0),'nullterminated'};
    
    headerspec.magic = char(zeros(1,8));
    headerspec.endian = zeros(1,4,'uint8');
    headerspec.filetype = int32(0);
    headerspec.title = tecplotstring;
    headerspec.numvar = int32(0);
    headerspec.varnames = {tecplotstring, 'specified','numvar'};
    
    zonespec.marker = single(0);
    zonespec.name = tecplotstring;
    zonespec.parent = int32(0);
    zonespec.id = int32(0);
    zonespec.solntime = double(0);
    zonespec.color = int32(0);
    zonespec.type = int16(0);
    zonespec.packing = int16(0);
    zonespec.isvarlocation = int32(0);
    zonespec.varlocation = {int32(0), 'specified','isvarlocation'};
    zonespec.isfaceneighbor = int32(0);
    zonespec.userconnections = int32(0);
    %only works after here if it's an ordered zone...
    zonespec.size = zeros(1,3,'int32');
    zonespec.isauxdata = int32(0);
    
    auxdataspec.name = tecplotstring;
    auxdataspec.format = int32(0);
    auxdataspec.value = tecplotstring;
    auxdataspec.ismoreauxdata = int32(0);
    
    fid = fopen(fn,'r');
    
    header = readbinaryrec(fid,headerspec,1);
    if (any(header.magic(1:5) ~= int32('#!TDV')))
        error('Not a TecPlot file.');
    end;
    
    zonespec.varlocation = {zeros(1,header.numvar,'int32'), 'specified','isvarlocation'};
    
    done = false;
    while (~done),
        %read zones
        zone1 = readbinaryrec(fid,zonespec,1);
        if (zone1.marker ~= 299.0)
            error('Cannot read file.  Structure may be corrupted.');
        end;
        
        done = zone1.isauxdata == 0;
        auxdata = auxdataspec;
        nauxdata = 0;
        while (~done)
            auxdata1 = readbinaryrec(fid,auxdataspec,1);
            done = auxdata1.ismoreauxdata == 0;
            nauxdata = nauxdata + 1;
            auxdata(nauxdata) = auxdata1;
        end;
        auxdata = auxdata(1:nauxdata);
        
        done = true;
    end;
    
    fclose(fid);
else
    fid = fopen(fn,'r');
    
    header.title = fgetregexptok(fid, 'TITLE\s+=\s+"([^"]*)"');
        
    ln = fgetl(fid);
    tok = regexp(ln, 'VARIABLES\s+=\s+"([^"]*)"', 'tokens', 'once');
    header.varnames = tok;
    done = false;
    while (~done)
        ln = fgetl(fid);
        tok = regexp(ln, '^\s*"([^"]*)"$', 'tokens', 'once');
        
        if (isempty(tok))
            done = true;
        else
            header.varnames = [header.varnames tok];
        end;
    end;
    nvars = length(header.varnames);

    [zone1.name] = fgetregexptok(fid, 'ZONE T="([^"]*)"', ln);
    
    fgetl(fid);
    
    M = fscanf(fid, ' I=%d, J=%d, K=%d, ZONETYPE=%s\n');
    zone1.size = M(1:3)';
    
    fgetl(fid);
    
    zone1.iscenter = false(1,nvars);
    ln = fgetl(fid);
    tok = regexp(ln, '(\[[0-9,-]+\])=CELLCENTERED', 'once','tokens');
    if (~isempty(tok))
        tok{1}(tok{1} == '-') = ':';
        
        indcenter = eval(tok{1});
        
        zone1.iscenter(indcenter) = true;
    end;
    
    ln = fgetl(fid);
    tok = regexp(ln, '([a-zA-Z]+)\s', 'tokens');
    zone1.datatype = tok;
    
    for i = 1:nvars,
        if (zone1.iscenter(i))
            sz = zone1.size - 1;
            sz(sz == 0) = 1;
        else
            sz = zone1.size;
        end;
        
        data1 = fscanf(fid, '%f ', prod(sz));
        
        data.(header.varnames{i}) = reshape(data1, sz);
    end;
end;
