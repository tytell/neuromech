function varargout = importsamrai(pathname,varargin)
% V = importsamrai(pathname)
%  or  [x,y,Vinterp,V] = importsamrai(pathname,'interpolaten',[100 50])
% Imports samrai data from files in the directory pathname.  Files should
% be nammed 'summary.samrai' and 'processor_clusterXXXXX.samrai', where the
% Xs are digits.  Also can interpolate onto an evenly spaced grid with the
% 'interpolaten' option, which specifies the number of interpolation
% positions in the x and y directions, respectively.
% 
% The output V is a structure containing one element for each grid region.
%   V.level_number: The grid resolution, from 1 (coarsest) usually to 4
%     (finest)
%   V.xlo, V.xup: Lower and upper bounds in x and y of the region
%   V.rows, V.cols: Number of rows and columns of the grid
% plus the data themselves, which are usually
%   V.U_0, V.U_1: u and v components of velocity
%   V.P: pressure
%   V.F0, V.F1: Artificial force field values from the IB points
%   V.Omega: vorticity
% everything in CGS units.

opt.debug = false;
opt.patchcol = 'kbrgcmy';
opt.ndigits = 5;
opt.interpolaten = [];
opt.interpolategrid = {};
opt.vars = {};
opt = parsevarargin(opt,varargin,2);

summaryfile = fullfile(pathname,'summary.samrai');
subfiles = getfilenames(fullfile(pathname,'processor_cluster*.samrai'));

%get the patch info
patchext = h5compound2struct(hdf5read(summaryfile, '/extents/patch_extents'));
patchmap = h5compound2struct(hdf5read(summaryfile, '/extents/patch_map'));

npatch = length(patchext);

%and the variable info
varnames0 = hdf5read(summaryfile,'/BASIC_INFO/var_names');
dx = hdf5read(summaryfile,'/BASIC_INFO/dx')';

varncomponents = hdf5read(summaryfile,'/BASIC_INFO/var_number_components');

%at the moment, we only import variables with a single component
good = varncomponents == 1;
varnames0 = varnames0(good);

%set up the structure
if (isempty(opt.vars)),
    nvar = length(varnames0);
    varnames = cell(size(varnames0));
    fieldnm = cell(size(varnames0));
    for i = 1:nvar,
        varnames{i} = varnames0(i).Data;
        tok = regexp(varnames0(i).Data,'::','split');
        fieldnm{i} = tok{end};
    end;
else
    fieldnm = opt.vars;
    
    shortvarnames = cell(size(varnames0));
    for i = 1:length(varnames0),
        tok = regexp(varnames0(i).Data,'::','split');
        shortvarnames{i} = tok{end};
    end;
    
    varnames = cell(size(opt.vars));
    nvar = length(opt.vars);
    for i = 1:nvar,
        k = find(strcmp(opt.vars{i},shortvarnames));
        
        if (length(k) == 1),
            varnames{i} = varnames0(k).Data;
        elseif (isempty(k)),
            % U_0 is the same as U_x in newer versions
            if strcmp(opt.vars{i},'U_0')
                k = find(strcmp('U_x', shortvarnames));
                if (length(k) == 1)
                    varnames{i} = varnames0(k).Data;
                end
            elseif strcmp(opt.vars{i},'U_1')
                k = find(strcmp('U_y', shortvarnames));
                if (length(k) == 1)
                    varnames{i} = varnames0(k).Data;
                end
            end
            
            if isempty(k)
                error('Could not find variable %s',opt.vars{i});
            end
        end;
    end;
end;

V = cell2struct(cell(npatch,nvar),fieldnm,2);

for f = 1:length(subfiles),
    %import the data from each patch
    fileid = H5F.open(subfiles{f},'H5F_ACC_RDONLY','H5P_DEFAULT');

    patchnum = find(cat(1,patchmap.processor_number) == f-1)';
    
    for i = patchnum,
        procname = sprintf('processor.%0*d',opt.ndigits,patchmap(i).processor_number);
        levelname = sprintf('level.%0*d',opt.ndigits,patchmap(i).level_number);
        patchname = sprintf('patch.%0*d',opt.ndigits,patchmap(i).patch_number);
        
        n = patchext(i).upper(1) - patchext(i).lower(1) + 1;
        m = patchext(i).upper(2) - patchext(i).lower(2) + 1;
        
        for j = 1:nvar,
            setname = ['/' procname '/' levelname '/' patchname '/' varnames{j}];
            
            [datasetid,err] = hdf5err(@H5D.open, fileid, setname);
            v = H5D.read(datasetid,'H5ML_DEFAULT','H5S_ALL','H5S_ALL','H5P_DEFAULT');
            
%             rng = [zeros(1,nd); makerow(datasetlen)];
%             rng = flatten(opt.range)';
% 
%             startind = rng(1:2:end);
%             endind = rng(2:2:end);
% 
%             len = endind - startind + 1;
% 
%             memspaceid = H5S.create_simple(nd,len,[]);
% 
%             %read the data
%             v = H5D.read(datasetid,'H5ML_DEFAULT',memspaceid,dataspaceid,...
%                 'H5P_DEFAULT');
%             
%             H5S.close(memspaceid);
            H5D.close(datasetid);

            v = reshape(v, [n m]);
            v = v';
            
            V(i).(fieldnm{j}) = v;
        end;
        V(i).level_number = patchmap(i).level_number;
        V(i).xlo = patchext(i).xlo;
        V(i).xup = patchext(i).xup;
        V(i).rows = m;
        V(i).cols = n;
    end;
    
    H5F.close(fileid);
end;

%sort the patches by level (not usually necessary)
[q,ord] = sort(cat(2,V.level_number));
V = V(ord);

%interpolate
if (~isempty(opt.interpolaten)),
    if (numel(opt.interpolaten) == 1),
        opt.interpolaten = opt.interpolaten([1 1]);
    end;
    [xi,yi,V1] = interpsamrai(V, 'numgrid',opt.interpolaten);

    
    varargout = {xi,yi,V1,V};
elseif (~isempty(opt.interpolategrid)),
    [xi,yi,V1] = interpsamrai(V, opt.interpolategrid{1},opt.interpolategrid{2});
    
    if (nargout == 4)
        varargout = {xi,yi,V1,V};
    else
        varargout = {V1,V};
    end;
else
    varargout = {V};
end;

if (opt.debug),
    clf;
    hold on;
    for i = 1:npatch,
        xlo = V(i).xlo;
        xup = V(i).xup;
        col = opt.patchcol(V(i).level_number+1);
        
        dx1 = dx(V(i).level_number+1,1);
        dy1 = dx(V(i).level_number+1,2);
        imagesc([xlo(1)+dx1/2 xup(1)-dx1/2], [xlo(2)+dy1/2 xup(2)-dy1/2], ...
            V(i).Omega);
        
        addplot([xlo(1) xlo(1) xup(1) xup(1) xlo(1)], ...
            [xlo(2) xup(2) xup(2) xlo(2) xlo(2)], col);
        drawnow;
    end;
    
    if (~isempty(opt.interpolaten)),
        addquiverc(xi,yi, V1.U_0, V1.U_1);
    end;
    
    hold off;
end;

function S = h5compound2struct(H)

fieldnm = H(1).MemberNames;

C = cell(length(H),length(fieldnm));
for i = 1:length(H),
    for j = 1:length(fieldnm),
        if (isa(H(i).Data{j},'hdf5.h5array'))
            C{i,j} = H(i).Data{j}.Data;
        else
            C{i,j} = H(i).Data{j};
        end;
    end;
end;
S = cell2struct(C,fieldnm,2);



