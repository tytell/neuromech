function batchanalyzesim(dirname, unattended)
% function batchanalyzesim(basedirname, unattended)
%   or     batchanalyzesim(datadirnames, unattended)
%   or     batchanalyzesim(filenames, unattended)
%
% Analyzes many data files.  Imports data from FORT.xx files, if necessary,
% then analyzes that data.  Parameters can be the base directory, which
% contains subdirectories that have the FORT.xx files, or a cell array of
% names of the data directories, or a cell array of existing Matlab data
% files.
%
% After the Fortran data is imported, it runs analyzesim on the file,
% producing another data file with the suffix '_analysis'.
%
% Produces a table containing useful results from each data file.

if (nargin == 1)
    unattended = false;
end;

inputyn('','clearsaved');
if (ischar(dirname) || (iscellstr(dirname) && exist(dirname{1},'dir'))),
    if (ischar(dirname)),
        datadirs = getdirectorynames(dirname);
    else
        datadirs = dirname;
    end;
    
    useexisting = unattended || inputyn('Use existing data files? ');
    
    basematfiles = cell(size(datadirs));
    samraibasedirs = cell(size(datadirs));
    goodfiles = true(size(datadirs));
    for i = 1:length(datadirs),
        [pathnm,name] = fileparts(datadirs{i});
        
        tok = regexpi(name,'(FORT)','split');
        name1 = strcat(tok{:});
        basematfiles{i} = fullfile(pathnm, ['sim' name1 '.mat']);
        samraibasedirs{i} = fullfile(pathnm, [tok{1} 'viz_IB2d']);
        if (~exist(samraibasedirs{i},'dir'))
            samraibasedirs{i} = fullfile(datadirs{i}, [tok{1} 'viz_IB2d']);
        end;
        
        if (useexisting),
            goodfiles(i) = exist(basematfiles{i},'file');
        else
            if (~exist(basematfiles{i},'file')),
                doimport = true;
                goodfiles(i) = false;
            elseif (~getvar('-file',basematfiles{i},'xn')),
                doimport = unattended || ...
                    inputyn(sprintf('Re-import data from directory %s (file seems old)? ', ...
                    datadirs{i}));
                goodfiles(i) = false;
            else
                if (unattended)
                    doimport = false;
                else
                    doimport = ~inputyn(sprintf('Use existing data in file %s (file seems OK)? ', ...
                        basematfiles{i}));
                end;
                goodfiles(i) = true;
            end;
            if (doimport),
                importsimdata(datadirs{i},basematfiles{i});
                goodfiles(i) = true;
            end;
        end;
    end;
elseif (iscellstr(dirname)),
    basematfiles = dirname;
    goodfiles = true(size(basematfiles));
    samraibasedirs = cell(size(basematfiles));
    for i = 1:length(basematfiles)
        [pathnm,name] = fileparts(basematfiles{i});
        
        samraibasedirs{i} = fullfile(pathnm, [name(4:end) 'viz_IB2d']);
    end;        
end;

figure(1);
clf;
set(gcf,'Color','w');

basematfiles = basematfiles(goodfiles);
samraibasedirs = samraibasedirs(goodfiles);
analysisfiles = cell(size(basematfiles));
rownames = cell(length(basematfiles),1);
tab = cell(length(basematfiles),12);
for i = 1:sum(goodfiles),
    fprintf('Getting constants for file %s...\n', basematfiles{i});
    
    [pn,fn] = fileparts(basematfiles{i});
    analysisfiles{i} = fullfile(pn,[fn '_analysis.mat']);
    
    freq = 1;
    viscosity = 0.01;
    sfo = 2.56e6;
    sfo2 = '';
    ps = 3;
    gridres = 32;
    isforcetaper = 'true';
    if (~getvar('-file',basematfiles{i}, 'freq','viscosity','sfo','sfo2','ps','gridres','isforcetaper',...
            '-keepundef') || ~inputyn('  Use existing constants from file?')),
        prompt = {'freq','viscosity','sfo','sfo2','ps','gridres','isforcetaper'};
        def = cellfun(@num2str,{freq,viscosity,sfo,sfo2,ps,gridres,isforcetaper},'UniformOutput',false);
        if (isempty(sfo2)),
            def{4} = '0.08*2.56e6';
        end;
        if (isforcetaper),
            def{7} = 'true';
        else
            def{7} = 'false';
        end;
        vals = inputdlg(prompt,'Simulation constants',1,def);
        
        vals = cellfun(@eval,vals,'UniformOutput',false);
        
        [freq,viscosity,sfo,sfo2,ps,gridres,isforcetaper] = vals{:};
        
        putvar('-file',basematfiles{i}, 'freq','viscosity','sfo','sfo2','ps','gridres','isforcetaper');
    end;
end;

for i = 1:sum(goodfiles),
    fprintf('Analyzing data in file %s...\n', basematfiles{i});

    goodfile = exist(analysisfiles{i},'file') &&  ...
        getvar('-file', analysisfiles{i}, 't','comspeed','issteady','steadycycle',...
        's','ampcont','amp','wavespeed','actspeed',...
        'wavelen','actlen','sfo','sfo2','viscosity','ps','freq','width', ...
        'worktot','workpos','workneg','workposact','worknegact', ...
        'Force');
        
    if (~goodfile || (~unattended && ~inputyn('Use existing results? '))),
        analyzesim(basematfiles{i}, 'outfile',analysisfiles{i}, ...
            'samraibasedir',samraibasedirs{i}, ...
            'savepressure','savevorticity');
        getvar('-file', analysisfiles{i}, 't','comspeed','issteady','steadycycle',...
            's','ampcont','amp','wavespeed','actspeed',...
            'wavelen','actlen','sfo','sfo2','viscosity','ps','freq','width',...
            'worktot','workpos','workneg','workposact','worknegact', ...
            'Force');
    end;
    
    [~,fn] = fileparts(basematfiles{i});
    rownames{i} = fn;
    
    len = s(end,1);
    
    tab(i,1:5) = {freq,sfo,sfo2,ps,viscosity};
    tab{i,6} = max(width(:,1)) / len;
    if (width(end-10,1) < 0.9*max(width(:,1))),
        tab{i,7} = 'yes';
    else
        tab{i,7} = 'no';
    end;
    tab{i,8} = nanmean(comspeed(issteady))/len;
    tab{i,9} = {t,comspeed/len,'k-'};
    
    tab{i,10} = nanmean(wavespeed(steadycycle:end))/len;
    tab{i,11} = tab{i,8} / tab{i,10};
    
    tab{i,12} = nanmean(flatten(amp(end-20:end,steadycycle:end)))/len;
    env = max(abs(ampcont(:,issteady))/len,[],2);
    tab{i,13} = {s([1 1:end end 1])/len, [0; env; 0; 0]/len, 'k-'};
    
    tab{i,14} = nanmean(flatten(actlen(end-20:end,steadycycle:end)))/len;
    tab{i,15} = nanmean(flatten(wavelen(end-20:end,steadycycle:end)))/len;
    tab{i,16} = nanmean(nanmean(wavelen(end-20:end,steadycycle:end)) ./ ...
        nanmean(actlen(end-20:end,steadycycle:end)));
    
    tab{i,17} = {s/len, nanmean(sum(workpos(:,steadycycle:end),3),2), 'k-', ...
        s/len, nanmean(sum(worknegact(:,steadycycle:end),3),2), 'r-', ...
        s/len, nanmean(sum(workneg(:,steadycycle:end),3),2), 'k-', ...
        s/len, nanmean(sum(worktot(:,steadycycle:end),3),2), 'b-'};
    
    tab{i,18} = {[0; linspace(0,1,size(Force.thrustt,1))'; 1], ...
        [0; nanmean(sum(Force.thrustt(:,steadycycle:end,:)+Force.thrustn(:,steadycycle:end,:),3),2); 0], 'g-', ...
        [0; linspace(0,1,size(Force.thrustt,1))'; 1], ...
        [0; -nanmean(sum(Force.dragt(:,steadycycle:end,:)+Force.dragn(:,steadycycle:end,:),3),2); 0], 'r-'};
    
    putvar basematfiles i tab;
    
    showtable(tab(1:i,:),'colnames',{'freq','sfo','sfo2','ps','mu','width','taper', ...
        'comspeed', '', 'wavespeed','slip', 'amp','env', 'actlen','wavelen','wave/act', ...
        'work','thrust'}, ...
        'rownames', rownames(1:i), 'align','lm', 'format','%.3g', 'colwidthmode','tight', ...
        'colgap',5,'rowgap',3, 'rowheight',20, 'xlim','linkcol','ylim','linkcol');
    drawnow;
end;

showtable(tab(1:i,:),'colnames',{'freq','sfo','sfo2','ps','mu','width','taper', ...
    'comspeed', '', 'wavespeed','slip', 'amp','env', 'actlen','wavelen','wave/act', ...
    'work','thrust'}, ...
    'rownames', rownames(1:i), 'align','lm', 'format','%.3g', 'colwidthmode','tight', ...
    'colgap',5,'rowgap',3, 'rowheight',20, 'xlim','linkcol','ylim','linkcol');

