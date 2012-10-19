function importsimdata(datadir,outfile, varargin)
%IMPORTSIMDATA  Imports Fortran output data from CFD simulations
%     importsimdata(datadir,outfile, opts...)
%
% datadir is the path to the directory containing the Fortran output files
% outfile is the name of a Matlab file for output
% opts include the names of all the different files

opt.muscforcefile = 'fort.35';
opt.totalforcefile = 'fort.33';
opt.xyfile = 'fort.31';
opt.uvfile = 'fort.32';
opt.actfile = 'fort.45';
opt.odeforcefile = 'fort.37';
opt.springbinfiles = {'SPRING.01','SPRING.02','SPRING.03','SPRING.04',...
    'SPRING.05','SPRING.06','SPRING.07','SPRING.08','SPRING.09','SPRING.10'};
%opt.springfiles = {'fort.101','fort.102','fort.103','fort.104',...
%    'fort.105','fort.106','fort.107','fort.108','fort.109','fort.110'};
opt.npt = 321;
opt.nlines = 4;

opt.springrectype = struct('i1',uint8(0),'j1',uint8(0),'i2',uint8(0),'j2',uint8(0),...
    'fx',double(0),'fy',double(0),'energy',double(0));

opt = parsevarargin(opt,varargin,2);

datafiles = {};

%xy coordinates of the IB points
fn = fullfile(datadir,opt.xyfile);
if (exist(fn,'file')),
    fprintf('Importing xy coordinates...\n');
    XY = dlmread(fn);
    datafiles = [datafiles fn];
    
    %check number of frames
    nfr = size(XY,1) / (opt.nlines * opt.npt);
    if (nfr ~= floor(nfr)),
        warning('importsimdata:nonintegerframes',...
                'Number of frames is non-integer in xy points file.  Truncating');
        ntrunc = floor(nfr) * opt.nlines * opt.npt;
        XY = XY(1:ntrunc,:);
        nfr = floor(nfr);
    end;
    fprintf('  %d frames.\n', nfr);
    
    nfr1 = nfr;
    xx = reshape(XY(:,1),[opt.npt opt.nlines nfr1]);
    yy = reshape(XY(:,2),[opt.npt opt.nlines nfr1]);
    
    out.xl = squeeze(xx(:,3,:));
    out.yl = squeeze(yy(:,3,:));
    out.xr = squeeze(xx(:,4,:));
    out.yr = squeeze(yy(:,4,:));
    out.xm = squeeze(xx(:,1,:));
    out.ym = squeeze(yy(:,1,:));
    
    out.xn = squeeze(xx(:,2,:));
    out.yn = squeeze(yy(:,2,:));
    
    %don't remove the ghost points
    %out.xl(end,:) = NaN;
    %out.yl(end,:) = NaN;
    %out.xr(end,:) = NaN;
    %out.yr(end,:) = NaN;
    %out.xn(end,:) = NaN;
    %out.yn(end,:) = NaN;
else
    warning('importsimdata:missingfile','No xy data file %s',fn);
end;

%uv velocities of the IB points
fn = fullfile(datadir,opt.uvfile);
if (~exist(fn,'file')),
    warning('importsimdata:missingfile','No uv data file %s',fn);
else    
    fprintf('Importing uv velocities...\n');
    UV = dlmread(fn);
    datafiles = [datafiles fn];
    
    %check number of frames
    nfr1 = size(UV,1) / (opt.nlines * opt.npt);
    if (nfr1 ~= floor(nfr1)),
        warning('importsimdata:nonintegerframes',...
                'Number of frames is non-integer in uv points file.  Truncating');
        ntrunc = floor(nfr1) * opt.nlines * opt.npt;
        UV = UV(1:ntrunc,:);
        nfr1 = floor(nfr1);
    end;
    if (nfr1 ~= nfr),
        fprintf('  ...but %d frames in uv velocities file...\n',nfr1);
    end;
    
    uu = reshape(UV(:,1),[opt.npt opt.nlines nfr1]);
    vv = reshape(UV(:,2),[opt.npt opt.nlines nfr1]);
    
    out.ul = squeeze(uu(:,3,:));
    out.vl = squeeze(vv(:,3,:));
    out.ur = squeeze(uu(:,4,:));
    out.vr = squeeze(vv(:,4,:));
    out.um = squeeze(uu(:,1,:));
    out.vm = squeeze(vv(:,1,:));

    out.un = squeeze(uu(:,2,:));
    out.vn = squeeze(vv(:,2,:));
    
    %don't remove the ghost points
    %out.ul(end,:) = NaN;
    %out.vl(end,:) = NaN;
    %out.ur(end,:) = NaN;
    %out.vr(end,:) = NaN;
    %out.un(end,:) = NaN;
    %out.vn(end,:) = NaN;
end;

%activation data
fn = fullfile(datadir,opt.actfile);
if (~exist(fn,'file')),
    warning('importsimdata:missingfile','No activation data file %s',fn);
else    
    fprintf('Importing activation data...\n');
    actdata = dlmread(fn);
    datafiles = [datafiles fn];
    
    %these are the indices for which activation is defined
    good = actdata(:,2) ~= 0;
    actdata = actdata(good,:);
    
    actpt = [min(actdata(:,2)) max(actdata(:,2))];
    nactpt = diff(actpt)+1;
    
    %check number of frames
    nfr1 = size(actdata,1) / nactpt;
    if (nfr1 ~= floor(nfr1)),
        warning('importsimdata:nonintegerframes',...
                'Number of frames is non-integer in activation file.  Truncating');
        ntrunc = floor(nfr1) * nactpt;
        actdata = actdata(1:ntrunc,:);
        nfr1 = floor(nfr1);
    end;
    if (nfr1 ~= nfr),
        fprintf('...but %d frames in activation file...\n',nfr1);
    end;
    
    %time - has values for each point and line, but they're the same, so only
    %save the values for each frame
    t = reshape(actdata(:,1),[nactpt nfr1]);
    out.t = t(1,:);
    
    out.actl = false(opt.npt,nfr1);
    out.actr = false(opt.npt,nfr1);
    
    out.actl(actpt(1):actpt(2),:) = reshape(actdata(:,3) > 0, [nactpt nfr1]);
    out.actr(actpt(1):actpt(2),:) = reshape(actdata(:,4) > 0, [nactpt nfr1]);
end;

%get the muscle force file
fn = fullfile(datadir,opt.muscforcefile);
if (~exist(fn,'file')),
    warning('importsimdata:missingfile','No muscle force data file %s',fn);
else    
    fprintf('Importing muscle force...\n');
    M = dlmread(fn);
    datafiles = [datafiles fn];
    
    nfr1 = size(M,1) / (opt.nlines * opt.npt);
    if (nfr1 ~= floor(nfr1)),
        warning('importsimdata:nonintegerframes',...
                'Number of frames is non-integer in muscle force file.  Truncating');
        ntrunc = floor(nfr1) * opt.nlines * opt.npt;
        M = M(1:ntrunc,:);
        nfr1 = floor(nfr1);
    end;
    if (nfr1 ~= nfr),
        fprintf('...but %d frames in muscle force file...\n',nfr1);
    end;
    
    %muscle force in x and y directions.  Only save the left and right sides,
    %because the force along the midline is zero
    fxmus = reshape(M(:,1),[opt.npt opt.nlines nfr1]);
    fymus = reshape(M(:,2),[opt.npt opt.nlines nfr1]);
    
    out.fxlmus = squeeze(fxmus(:,3,:));
    out.fylmus = squeeze(fymus(:,3,:));
    out.fxrmus = squeeze(fxmus(:,4,:));
    out.fyrmus = squeeze(fymus(:,4,:));
end;

%total force data
fn = fullfile(datadir,opt.totalforcefile);
if (~exist(fn,'file')),
    warning('importsimdata:missingfile','No total force data file %s',fn);
else    
    fprintf('Importing total force...\n');
    T = dlmread(fn);
    datafiles = [datafiles fn];
    
    %check number of frames
    nfr1 = size(T,1) / (opt.nlines * opt.npt);
    if (nfr1 ~= floor(nfr1)),
        warning('importsimdata:nonintegerframes',...
                'Number of frames is non-integer in total force file.  Truncating');
        ntrunc = floor(nfr1) * opt.nlines * opt.npt;
        T = T(1:ntrunc,:);
        nfr1 = floor(nfr1);
    end;
    if (nfr1 ~= nfr),
        fprintf('...but %d frames in total force file...\n',nfr1);
    end;
    
    fxtot = reshape(T(:,1),[opt.npt opt.nlines nfr1]);
    fytot = reshape(T(:,2),[opt.npt opt.nlines nfr1]);
    out.fxltot = squeeze(fxtot(:,3,:));
    out.fyltot = squeeze(fytot(:,3,:));
    out.fxrtot = squeeze(fxtot(:,4,:));
    out.fyrtot = squeeze(fytot(:,4,:));
    
    out.fxmtot = squeeze(fxtot(:,1,:));
    out.fymtot = squeeze(fytot(:,1,:));
    out.fxntot = squeeze(fxtot(:,2,:));
    out.fyntot = squeeze(fytot(:,2,:));
    
    %don't remove the ghost points
    %out.fxltot(end,:) = NaN;
    %out.fyltot(end,:) = NaN;
    %out.fxrtot(end,:) = NaN;
    %out.fyrtot(end,:) = NaN;
    %out.fxntot(end,:) = NaN;
    %out.fyntot(end,:) = NaN;
end;

%ode solution data
fn = fullfile(datadir,opt.odeforcefile);
if (~exist(fn,'file')),
    warning('importsimdata:missingfile','No ODE force data file %s',fn);
else    
    fprintf('Importing ODE solution data...\n');
    odedata = dlmread(fn);
    datafiles = [datafiles fn];
    
    %these are the indices for which activation is defined
    if ((size(odedata,2) == 2) || ...
            ((size(odedata,2) == 3) && all(odedata(:,3) == 0))),
        odept = [min(odedata(:,1)) max(odedata(:,1))];
        nodept = diff(odept)+1;
        
        nside = 2;
        nfr1 = size(odedata,1) / (2 * nodept);

        if (nfr1 ~= floor(nfr1)),
            warning('importsimdata:nonintegerframes',...
                    'Number of frames is non-integer in ODE force file.  Truncating');
            ntrunc = floor(nfr1) * 2 * nodept;
            odedata = odedata(1:ntrunc,:);
            nfr1 = floor(nfr1);
        end;
        if (nfr1 ~= nfr),
            fprintf('...but %d frames in ODE force file...\n',nfr1);
        end;
        
        fmusP = zeros(opt.npt,nside,nfr1);
        fmusP(odept(1):odept(2),:,:) = reshape(odedata(:,2), [nodept nside nfr1]);
    else
        odept = [min(odedata(:,2)) max(odedata(:,2))];
        nodept = diff(odept)+1;
        
        %activation side
        side = odedata(:,1);
        nside = max(side);
        
        %check number of frames
        nfr1 = size(odedata,1) / (nside * nodept);
        if (nfr1 ~= floor(nfr1)),
            error('Number of frames is non-integer in activation data file');
        elseif (nfr1 ~= nfr),
            fprintf('...but %d frames in activation file...\n',nfr1);
        end;
        
        fmusP = zeros(opt.npt,nside,nfr1);
        fmusP(odept(1):odept(2),:,:) = reshape(odedata(:,3), [nodept nside nfr1]);
    end;
    
    out.fmusPl = squeeze(fmusP(:,1,:));
    out.fmusPr = squeeze(fmusP(:,2,:));
end;

%total force data
if (~isempty(opt.springbinfiles) && ~exist(fullfile(datadir,opt.springbinfiles{1}),'file')),
    disp('No spring force data.  Skipping...');
else    
    fprintf('Importing spring force...\n');
    
    nspring = length(opt.springbinfiles);
    for i = 1:length(opt.springbinfiles),
        fn = fullfile(datadir,opt.springbinfiles{i});
        datafiles = [datafiles fn];
        
        S1 = readfortranrecs(fn,opt.springrectype);
        %stupid error - we wrote out the points from 1:320ish in an 8bit
        %field, so it wraps around at 255.  Look for exact jumps down by
        %255
        off = NaN(size(S1.i1));
        iswrap = [false; (S1.i1(1:end-1) == 255) & (S1.i1(2:end) == 0)];
        off(iswrap) = 256;
        isreset = [false; (S1.i1(1:end-1) > S1.i1(2:end))];
        isreset = isreset & ~iswrap;
        off(isreset) = 0;
        off(1) = 0;
        off = fillmat(off);
        S1.i1 = uint16(S1.i1) + uint16(off);
        
        off = NaN(size(S1.i2));
        iswrap = [false; (S1.i2(1:end-1) == 255) & (S1.i2(2:end) == 0)];
        off(iswrap) = 256;
        isreset = [false; (S1.i2(1:end-1) > S1.i2(2:end))];
        isreset = isreset & ~iswrap;
        off(isreset) = 0;
        off(1) = 0;
        off = fillmat(off);
        S1.i2 = uint16(S1.i2) + uint16(off);
        
        nspringpt = range(S1.i1)+1;
        nfr1 = length(S1.i1) / double(nspringpt);
        if (nfr1 ~= floor(nfr1)),
            warning('importsimdata:nonintegerframes',...
                'Number of frames is non-integer in spring force file %s.  Truncating', ...
                opt.springbinfiles{i});
            ntrunc = floor(nfr1) * nspringpt;

            fieldnms = fieldnames(S1);
            for j = 1:length(fieldnms),
                fnm = fieldnms{j};
                S1.(fnm) = S1.(fnm)(1:ntrunc);
            end;
            nfr1 = floor(nfr1);
        end;
        if (nfr1 ~= nfr),
            fprintf('...but %d frames in total force file...\n',nfr1);
        end;

        fieldnms = fieldnames(S1);
        for j = 1:length(fieldnms),
            fnm = fieldnms{j};
            S1.(fnm) = reshape(S1.(fnm),[nspringpt nfr1]);
        end;
        
        if (i == 1),
            out.fxsp = zeros(opt.npt,nfr1,opt.nlines,nspring,2);
            out.fysp = zeros(opt.npt,nfr1,opt.nlines,nspring,2);
            out.energy = zeros(opt.npt,nfr1,nspring);
        end;
        
        for k = 1:nspringpt,
            out.fxsp(S1.i1(k,1),:,S1.j1(k,1),i,1) = S1.fx(k,:);
            out.fxsp(S1.i2(k,1),:,S1.j2(k,1),i,2) = -S1.fx(k,:);
            out.fysp(S1.i1(k,1),:,S1.j1(k,1),i,1) = S1.fy(k,:);
            out.fysp(S1.i2(k,1),:,S1.j2(k,1),i,2) = -S1.fy(k,:);
        end;
        
        out.energy(S1.i1(:,1),:,i) = S1.energy;
    end;
end;

issimname = false;
if (~strcmp(datadir,'.')),
    [~,dirname] = fileparts(datadir);
    
    tok = regexp(dirname, '-', 'split');
    if (length(tok) == 4),
        out.SIMdate = str2double(tok(1:3));

        numstr = regexp(tok{4}, '\d{4,}', 'match', 'once');
        if (~isempty(numstr)),
            out.SIMnum = str2double(numstr);
        end;
        out.SIMname = tok{4};
    
        id = ['SIM' dirname];
        id(id == '-') = '_';
        id = genvarname(id);
        out.(id) = 1;
        issimname = true;
    end;
end;
if (~issimname),
    [~,fn] = fileparts(outfile);
    id = upper(fn);
    id(id == '-') = '_';
    id = genvarname(id);
    out.(id) = 1;
    issimname = true;
end;

fields = fieldnames(out);
nfr = structfun(@(x) (size(x,2)), out);
ufr = unique(nfr);

if ((length(ufr) >= 2) && (ufr(end-1) ~= 1)),
    if (ufr(end-1) == ufr(end)-1),
        k = 2:ufr(end);
    else
        %if anything has one extra frame, it means that there's output at
        %time 0.  we'll discard that
        
        fac = floor(ufr(end) / ufr(end-1));
        
        fprintf('Time step in some variables seems to be %dx in others\n', fac);
        
        k = 1:fac:ufr(end);
        k = k(1:ufr(end-1));
        out.fac = fac;
    end;
    
    for i = 1:length(fields),
        if (size(out.(fields{i}),2) == ufr(end)),
            v1 = out.(fields{i});
            nm0 = [fields{i} '0'];
            out.(nm0) = out.(fields{i});

            sz = size(v1);
            v1 = permute(v1,[2 1 3:ndims(v1)]);
            v1 = flatten(v1,2:ndims(v1));
            v1 = v1(k,:);
            v1 = reshape(v1,[length(k) sz([1 3:end])]);
            v1 = ipermute(v1,[2 1 3:ndims(v1)]);
            out.(fields{i}) = v1;
        end;
    end;
end;

out.HGREV = savehgrev([],'datafile',datafiles);

save(outfile,'-struct','out');



