function [piv,names,units] = analyzePIVload(file)

piv = [];
names = [];
units = [];

if ((nargin == 0) || ~exist(file,'file')),
    filter = {'*.vec','Insight files (*.vec)'; ...
        '*.dat;*.vec','TecPlot files (*.dat,*.vec)'; ...
              '*.imx;*.vec','DaVis 6 files (*.imx,*.vec)'; ...
              '*.im7;*.vc7','DaVis 7 files (*.im7,*.vc7)'; ...
              '*.nc','NetCDF files (*.nc)'; ...
              '*.mat','Matlab data files (*.mat)'};
    filters = {'insight','tecplot','davis','davis','netcdf','matlab'};

    [file,path,filtind] = uigetfile(filter, 'Select File(s)', 'MultiSelect','on');
    if filtind == 0
        return;
    else
        where = filters{filtind};
    end
else
    [~,~,ext] = fileparts(file);
    switch lower(ext)
        case '.dat'
            where = 'tecplot';
        case {'.imx','.vec','.im7','.vc7'}
            where = 'davis';
        case '.nc'
            where = 'netcdf';
        case '.mat'
            where = 'matlab';
    end
end

if (~iscell(file))
    file = fullfile(path,file);
    [files,frames] = apMatchFileSeq(file);
    piv.frames = frames;
else
    files = file;
    frames = apGetFrameNumbers(files);
end

switch where,
    case 'matlab',
        errordlg('Matlab files currently not supported.');
        return;
        
    case 'insight',
        [x1,y1,u1,v1,~,~,insightdata] = loadInsight(files{1});
        piv.x = x1;
        piv.y = y1;
        piv.u = u1;
        piv.v = v1;
        piv.w = [];  % right now only load in u,v vectors
        sz = size(u1);
        
        N = length(files);
        if (N > 20),
            timedWaitBar(0,'Loading frames...');
        end;
        for i = 2:N,
            [x1,y1,u1,v1] = loadInsight(files{i});
            if (any(size(u1) ~= sz)),
                lasterr('Insight files have different numbers of vectors.');
                uiresume;
            end;
            if (any(size(x1) ~= size(u1))),
                [x1,y1] = meshgrid(x1,y1);
            end;
            piv.x(:,:,i) = x1;
            piv.y(:,:,i) = y1;
            piv.u(:,:,i) = u1;
            piv.v(:,:,i) = v1;
            
            if ((N > 20) && ~timedWaitBar(i/length(files))),
                break;
            end;
        end;
        if (N > 20),
            timedWaitBar(1);
        end;
        if (i ~= length(files)),
            warndlg(sprintf('Only loaded %d frames of %d.',i,N));
        end;
        piv.files = files;
        piv.frames = frames;
        names = insightdata.VARIABLES;
        
        units.x = insightdata.LengthUnit;
        units.vel = {insightdata.LengthUnit, insightdata.TimeUnit};
        
        data.vnames = {'Insight'};
        data.vsizes = size(piv.u)';
        data.piv = piv;
        data.pivunits = units;
        w1 = [];
        
    case 'davis',
        [x1,y1,u1,v1,units] = readimxvec(files{1});
        if (any(size(x1) ~= size(u1))),
            [x1,y1] = meshgrid(x1,y1);
        end;
        piv.x = x1;
        piv.y = y1;
        piv.u = u1;
        piv.v = v1;
        piv.w = [];  % right now only load in u,v vectors
        sz = size(u1);
        
        N = length(files);
        if (N > 20),
            timedWaitBar(0,'Loading frames...');
        end;
        for i = 2:N,
            [x1,y1,u1,v1,units] = readimxvec(files{i});
            if (any(size(u1) ~= sz)),
                lasterr('IMX files have different numbers of vectors.');
                uiresume;
            end;
            if (any(size(x1) ~= size(u1))),
                [x1,y1] = meshgrid(x1,y1);
            end;
            piv.x(:,:,i) = x1;
            piv.y(:,:,i) = y1;
            piv.u(:,:,i) = u1;
            piv.v(:,:,i) = v1;
            
            if ((N > 20) && ~timedWaitBar(i/length(files))),
                break;
            end;
        end;
        if (N > 20),
            timedWaitBar(1);
        end;
        if (i ~= length(files)),
            warndlg(sprintf('Only loaded %d frames of %d.',i,N));
        end;
        piv.files = files;
        piv.frames = frames;
        
        data.vnames = {'DaVis'};
        data.vsizes = size(piv.u)';
        data.piv = piv;
        data.pivunits = units;
        
        [k,q,tok] = regexp(units.x,'\[(.+)\]');
        if (~isempty(k)),
            units.pos = units.x(tok{1}(1,1):tok{1}(1,2));
        end;
        [k,q,tok] = regexp(units.vel,'\[?(.+)/(.+)\]?');
        if (~isempty(k)),
            units.vel{1} = units.vel(tok{1}(1,1):tok{1}(1,2));
            units.vel{2} = units.vel(tok{1}(2,1):tok{1}(2,2));
        end;
        
    case 'tecplot',
        i = 1;
        N = length(files);
        [x1,y1,u1,v1,w1,names] = apReadDat(files{i});
        piv.x = zeros(size(x1,1), size(x1,2), N);
        piv.y = zeros(size(y1,1), size(y1,2) ,N);
        piv.u = zeros(size(u1,1), size(u1,2), N);
        piv.v = zeros(size(v1,1), size(v1,2), N);
        if ~isempty(w1),
            piv.w = zeros(size(w1,1), size(w1,2), N);
        else
            piv.w = [];
        end
        
        piv.x(:,:,i) = x1;
        piv.y(:,:,i) = y1;
        piv.u(:,:,i) = u1;
        piv.v(:,:,i) = v1;
        if ~isempty(w1),
            piv.w(:,:,i) = w1;
        end
        
        timedWaitBar(0,'Loading...');
        for i = 2:N,
            [x1,y1,u1,v1,w1,names] = apReadDat(files{i});
            piv.x(:,:,i) = x1;
            piv.y(:,:,i) = y1;
            piv.u(:,:,i) = u1;
            piv.v(:,:,i) = v1;
            if ~isempty(w1),
                piv.w(:,:,i) = w1;
            end  % pm
            
            if (~timedWaitBar(i/N)),
                piv = [];
                names = [];
                units = [];
                return;  % cancel the load
            end;
        end;
        timedWaitBar(1);  % close timedWaitBar
        k = find((piv.u == 0) & (piv.v == 0));
        piv.u(k) = NaN;
        piv.v(k) = NaN;
        if ~isempty(w1),
            piv.w(k) = NaN;
        end
        
    case 'netcdf',
        errordlg('NetCDF not supported yet.');
        return;
end;

piv.files = files;
piv.fileType = where;
if ~isempty(w1),
    piv.isStereo = true;  % 3 components of velocity (stereo piv)
else
    piv.isStereo = false; % only planar velocity loaded
end


% ------------------------------------------------------------
function [x,y,u,v,w,vnames] = apReadDat(file)

x = [];
y = [];
u = [];
v = [];
w = [];

fid = fopen(file);
ln = fgetl(fid);                        % discard TITLE line

ln = fgetl(fid);                        % VARIABLES line
if (isempty(strfind(ln,'VARIABLES'))),
    fclose(fid);
    return;
end;
[a,b] = regexpi(ln,'"(\w+)"');
nvar = length(a);
for i = 1:nvar,
    vnames{i} = ln(a(i)+1:b(i)-1);
end;

ln = fgetl(fid);
if (isempty(strfind(ln,'ZONE'))),
    fclose(fid);
    return;
end;
[a,b,tok] = regexpi(ln,'(\w+)=(\w+)');
for i = 1:length(a),
    switch ln(tok{i}(1,1):tok{i}(1,2)),
        case 'I',
            r = str2num(ln(tok{i}(2,1):tok{i}(2,2)));
        case 'J',
            c = str2num(ln(tok{i}(2,1):tok{i}(2,2)));
    end;
end;

if (nvar == 4),  % x,y,u,v
    data = fscanf(fid,'%f %f %f %f\n',[4 Inf]);
else % nvar == 5  x,y,u,v,w
    data = fscanf(fid,'%f %f %f %f %f\n',[5 Inf]);
end;

fclose(fid);

data = data';
x = reshape(data(:,1),[r c])';
y = reshape(data(:,2),[r c])';
u = reshape(data(:,3),[r c])';
v = reshape(data(:,4),[r c])';
if (nvar == 5),
    w = reshape(data(:,5),[r c])';
end;

% ------------------------------------------------------------
function [files,frames] = apMatchFileSeq(file)

num1 = [];

tok = regexpi(file,'(.*[a-z])([0-9]+)(\.\w{1,3})$','tokens','once');
if (~isempty(tok)),
    base = tok{1};
    num1 = str2double(tok{2});
    ext = tok{3};
    expr = '.*\D([0-9]+)\.w{1,3}';
else
    %insight format
    tok = regexpi(file,'^(.*\D)([0-9]+)(\.T\d+\.\D\d+\.\P\d+\.\H\d+\.[LR]\.\w{1,3})$','tokens','once');
    
    if (~isempty(tok))
        expr = '.*\D([0-9]+)\.T';
        base = tok{1};
        num1 = str2double(tok{2});
        ext = tok{3};
    end
end

if (~isempty(num1))
    names = getfilenames(strcat(base,'*',ext));
    tok = regexp(names, expr, 'tokens','once');
    fnum = zeros(size(tok));
    for i = 1:length(tok)
        fnum(i) = str2double(tok{i});
    end
    [fnum,ord] = sort(fnum);
    names = names(ord);
    
    done = 0;
    num1 = num2str(num1);
    num2 = num2str(fnum(end));
    answer = {num1, num2};
    while (~done),
        answer = inputdlg({'Start:','End:'},'Select file numbers', ...
            1, answer);
        
        if (~isempty(answer)),
            i1 = find(fnum == str2double(answer{1}));
            i2 = find(fnum == str2double(answer{2}));
            if (isempty(i1) || isempty(i2)),
                waitfor(errordlg('Incorrect file number'));
            else
                done = 1;
            end;
        else
            done = 1;
        end;
    end;
    
    if (~isempty(answer)),
        files = names(i1:i2);
    end;
    
    frames = i1:i2;
else
    files = {file};
    frames = 1;
end;

