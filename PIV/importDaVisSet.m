function importDaVisSet(outfile,infiles)

if (nargin < 2),
    [fn,pn] = uigetfile({'*.vec;*.vc7','DaVis vector files (*.vec,*.vc7)'; ...
        '*.*', 'All files'}, 'Choose DaVis input files', ...
        'MultiSelect','on');
    
    %attach the path name to every file
    for i = 1:length(fn),
        infiles{i} = fullfile(pn,fn{i});
    end;
    
    if (nargin == 0),
        [fn,pn] = uiputfile({'*.mat','Matlab files (*.mat)'; ...
            '*.*', 'All files'}, 'Choose Matlab output file');
        
        outfile = fullfile(pn,fn);
    end;
end;

nfiles = length(infiles);

%open the first file so that we know how big everything should be
A = readimx(infiles{1});

nx = A.Nx;
ny = A.Ny;

rx = 1:nx;
ry = 1:ny;

x = (rx-1)*A.Grid*A.ScaleX(1) + A.ScaleX(2);
y = (ry-1)*A.Grid*A.ScaleY(1) + A.ScaleY(2);

filetype = A.IType;
switch A.IType,
    case 0,
        error('Cannot read image files.');
        
    case {1,3,5},
        if (A.IType == 5),
            nvec = 3;                   % 3 component vectors
        else
            nvec = 2;                   % 2 component vectors
        end;
        
    otherwise,
        error('Unsupported DaVis file type (#%d)',A.IType);
end;

doublesize = 8;                         % number of bytes in a double
totalsize = nx * ny * nfiles * nvec * doublesize;

if (totalsize > 10*2^20),                % greater than 10Mb?
    warning('Total file size will be %d Mb',round(totalsize/(2^20)));
end;

%initialize everything
u = zeros(ny,nx,nfiles);
v = zeros(ny,nx,nfiles);
if (nvec == 3),
    w = zeros(ny,nx,nfiles);
end;

for fnum = 1:nfiles,
    A = readimx(infiles{fnum});
    
    npages = size(A.Data,2)/ny;
    %error if we can't divide A.Data into a discrete number of pages
    if (mod(npages,1) ~= 0),
        error('Non-integer number of pages in Data structure.');
    end;

    data = double(reshape(A.Data,[nx ny npages]));

    %first page in data is a 1-based index into successive pages for which
    %vector DaVis chose.  Last pages are postprocessed vectors (given indices
    %4 or 5).  But since there are a max of 4*nvec pages, we have to
    %truncate at pageind == 3.  Don't ask me why
    pageind = data(:,:,1);
    pageind(pageind >= 4) = 4;
    %remember where the empty vectors are
    bad = pageind == 0;
    %but set the index so it won't give an error
    pageind(bad) = 1;

    %add 2 to skip over the pageind data itself
    pageind = (pageind-1)*nvec + 2;

    [i,j] = ndgrid(1:nx,1:ny);
    ind = sub2ind(size(data),i,j,pageind);
    u1 = data(ind);

    ind = sub2ind(size(data),i,j,pageind+1);
    v1 = data(ind);

    u1(bad) = NaN;
    v1(bad) = NaN;

    if (nvec == 3),
        ind = sub2ind(size(data),i,j,pageind+2);
        w1 = data(ind);
        w1(bad) = NaN;
    end;
    
    %flip and transpose from DaVis's normal layout until we end up with the
    %usual order for Matlab
    u1 = u1'*A.ScaleI(1) + A.ScaleI(2);
    v1 = v1'*A.ScaleI(1) + A.ScaleI(2);
    u1 = flipud(u1);
    v1 = flipud(v1);

    u(:,:,fnum) = u1;
    v(:,:,fnum) = v1;
    
    if (nvec == 3),
        w1 = w1'*A.ScaleI(1) + A.ScaleI(2);
        w1 = flipud(w1);
        
        w(:,:,fnum) = w1;
    end;
end;

%flip y and make a normal meshgrid type position matrix
y = y(end:-1:1);
[x,y] = meshgrid(x,y);

%check to see if we have the usual DaVis situation with velocities in m/s
%and positions in mm
%if so, convert velocities to mm/s so that everything is consistent
if (((regexp(A.UnitX,'^\[?mm\]?$') & ...
        regexp(A.UnitI,'^\[?m/s\]?$'))) | ...
        ((regexp(A.LabelX,'^\[?mm\]?$') & ...
        regexp(A.LabelI,'^\[?m/s\]?$')))),
    u = u*1000;
    v = v*1000;
    if (nvec == 3),
        w = w*1000;
    end;
end;

if (nvec == 2),
    savevars = {'x','y','u','v','infiles'};
else
    savevars = {'x','y','u','v','w','infiles'};
end;

save(outfile,savevars{:});
