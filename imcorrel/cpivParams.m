function [imageData,outname,gridData, optstruct] = cpivParams(params)

[imageData, params] = cpivProcessImageParams(params);
if (ischar(params{1})),
    outname = params{1};
    params = params(2:end);
    
    fid = fopen(outname,'r');
    if (fid ~= -1),
        fclose(fid);
        if (~inputyn(sprintf('File %s already exists.  Overwrite? ',outname),0)),
            fprintf('Cannot overwrite file.  Skipping.\n');
            optstruct = [];
            outname = [];
            gridData = [];
            return;
        end;
    end;
else
    outname = [];
end;

% find the region
[rgn,params] = getOption(params, 'Region', []);

% Get the grid
[gridData, params] = cpivProcessGridParams(params, imageData.Width, imageData.Height, rgn);

% Deal with object issues
[obj,params] = getOption(params, 'Object', []);
[overlap,params] = getOption(params, 'Overlap', 0.25);
[objinterp,params] = getOption(params, 'ObjectMode', 'none');
gridData.Object = obj;
gridData.Overlap = overlap;
gridData.ObjInterp = objinterp;

[interpmode,params] = getOption(params, 'InterpolationMode', '');
[snrmin,params] = getOption(params, 'HartCorrection', 1.0);
[savemode,params] = getOption(params, 'SaveMode', 1);
[uvbkgnd,params] = getOption(params, 'BackgroundVelocity', [0 0]);
[dosymm,params] = getOption(params, 'Symmetric', 0);

% set the default interpolation mode (agw if number of vectors > 1000, sts
% otherwise).  This is because sts gets *very* slow for large numbers of
% vectors
if (isempty(interpmode)),
    if (max(prod(gridData.NPt(1:end-1,:),2)) > 1000),
        interpmode = 'agw';
    else
        interpmode = 'sts';
    end;
end;
optstruct.InterpolationMode = interpmode;
if (length(snrmin) == 1),
    snrmin = repmat(snrmin,[gridData.NPasses 1]);
end;
optstruct.HartCorrection = snrmin;
optstruct.SaveMode = savemode;
optstruct.Symmetric = dosymm;
if (dosymm),
    warning('Symmetric option is very beta.  Use at your own risk.');
    
    % have to have an odd number of rows and columns in the correlation
    % matrix to do the symmetric processing
    gridData.SearchSize = gridData.SearchSize + ...
        1 - mod(gridData.SearchSize-gridData.WindowSize, 2);
end;

% set up transformation data for background flow
if (all(uvbkgnd == 0)),
    uvtform = [];
elseif (all(uvbkgnd-floor(uvbkgnd) == 0)),		% integer
    uvtform = -uvbkgnd;
else											% fraction
    uvtform = maketform('affine',[1 0 0; 0 1 0; -uvbkgnd 1]);
end;
optstruct.UVBkgnd = uvbkgnd;
optstruct.I2tform = uvtform;

if (~isempty(params)),
    error('Unrecognized parameters.');
end;

%------------------------------%
% utility function for scanning for options
function [val,options] = getOption(options, name,default)

q = find(cellfun('isclass',options, 'char'));
ind = strmatch(name,options(q));
ind = q(ind);

if (isempty(ind)),
    val = default;
elseif (ind + 1 <= length(options)),
    ind = ind(1);
    
    val = options{ind+1};
    options = {options{1:ind-1} options{ind+2:end}};
else
    error(sprintf('You must supply a value after ''%s''.', name));
end;
