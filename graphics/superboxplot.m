function superboxplot(varargin)
% function superboxplot(x,y,gp, opts...)
%
% Does a standard statistical box plot (like boxplot) but includes the x
% position.  x and y can be the same size, or x can be a column and y can be
% a matrix with the same number of rows as x.  gp must be the same size as y
% and indicates variable grouping beyond those that have the same x value.
% Boxes show the statistical distribution of points with the same x value
% and group value.  Boxes have a white line at the median and reach from 25
% to 75 quartile.  The notch indicates an approximate 95% confidence
% interval, and whiskers are 1.5 times the interquartile range long.  Any
% points that fall outside of the whisker are shown as individual points.
%
% Optional parameters:
%    'width' (or 'w') - Width of the boxes.  superboxplot attempts to set a reasonable
%      value for the width, but it doesn't always succeed.
%    'color' (or 'col') - Color order for boxes.  Uses the default color order if this
%      parameter isn't passed.
%    'fill' - false if you want the boxes not to be filled.
%    'xoffset' (or xoff or offx) - Offset data in different groups, but the
%      same x, by this amount.
%    'notch' (or 'n') - false if you don't want to show the notch that indicates
%      statistical significance.
%    'outliers' (or out) - false if you don't want to show outliers.

% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2011, Eric Tytell <tytell at jhu dot edu>

opt.width = [];
opt.fill = true;
opt.open = logical([]);
opt.color = '';
opt.xoffset = [];
opt.notch = true;
opt.outliers = true;
opt.marker = '.';

optsyn = { ...
    {'width',{'w'}}, ...
    {'fill',{'f'}}, ...
    {'open',{'o'}}, ...
    {'color',{'col'}}, ...
    {'xoffset',{'xoff','offx','offsetx'}}, ...
    {'notch',{'n'}}, ...
    {'outliers',{'out'}}};
    
[opt,arg] = parsevarargin(opt,varargin, 'leaveunknown','typecheck',false, ...
    'synonyms',optsyn);
    
if (length(arg) == 3)
    [x,y,gp] = deal(arg{:});
elseif (length(arg) == 2)
    [x,y] = deal(arg{:});
    gp = [];
end;

if (~isempty(opt.open) && opt.open),
    opt.fill = false;
end;

%look for groups
showgroups = true;
if ((ndims(y) > ndims(x)) || any(size(x) ~= size(y))),
    if (ndims(y) > 2),
        szx = size(x);
        szy = size(y);
        szx = [szx ones(1,length(szy)-length(szx))];
        
        m = prod(szx(szx == szy));
        n = prod(szy(szy ~= szx));
        
        x = reshape(x,[m 1]);
        y = reshape(y,[m n]);
    end;
    
    x = shiftdim(x);
    if (length(x) == size(y,2)),
        y = y';
    end;
    
    x = repmat(x,[1 size(y,2)]);
    if (isempty(gp)),
        gp = repmat(1:size(y,2),[size(x,1) 1]);
        showgroups = false;
    end;
elseif (all(size(x) == size(y)) && (size(y,1) == 1)),
    %deal with row vectors
    y = y';
    x = x';
end;
if (isempty(gp)),
    gp = repmat(1:size(y,2),[size(y,1) 1]);
    showgroups = false;
end;

%get rid of NaNs
good = isfinite(x) & isfinite(y);
x = x(good);
y = y(good);
[x,ind] = sort(x(:));
y = y(ind);

if (~isempty(gp)),
    gp = gp(good);
    gp = gp(ind);
    
    [gps,~,ind] = unique(gp(:));
    ngp = length(gps);
else
    gp = ones(size(y));
    gps = 1;
    ngp = 1;
end;

%look for categories in x
xind = find(diff(x) ~= 0)+1;
nx = length(xind)+1;
if (nx > 0.1*length(x)),
    warning('superboxplot:notcategorical','x variable does not seem to be categorical.');
end;
fprintf('%d x values detected: ',nx);
fprintf('%f ',x(xind-1));
fprintf('%f\n',x(end));

if (showgroups),
    fprintf('%d groups detected: ',length(gps));
    if (isnumeric(gps)),
        fprintf('%f ',gps);
    else
        fprintf('%c ',gps);
    end;
    fprintf('\n');
end;

%handle a single x value
rx = range(x(:));
if (rx == 0),
    rx = 1.1*ngp;
end;

%offset group values within a single x
if (isempty(opt.xoffset)),
    dx0 = rx/(1.1*nx*ngp);
    d = 1:ngp;
    d = d - (d(end)+d(1))/2;
    opt.xoffset = d*dx0;
elseif (numel(opt.xoffset) == 1),
    dx0 = opt.xoffset;
    d = 1:ngp;
    d = d - (d(end)+d(1))/2;
    opt.xoffset = d*dx0;
end;

%set the width of the boxes
if (isempty(opt.width)),
    if (opt.xoffset == 0),
        ntot = nx;
    else
        ntot = nx*ngp;
    end;
    opt.width = rx/(1.2*ntot);
end;

%color
if (isempty(opt.color)),
    opt.color = get(gca,'ColorOrder');
end;
if (size(opt.color,1) < length(gps)),
    opt.color = repmat(opt.color,...
                           [ceil(length(gps)/size(opt.color,1)) 1]);
end;

%handle "hold on" settings
switch get(gca,'NextPlot'),
 case 'replace',
  cla reset;
  unhold = true;
 case 'replacechildren',
  cla;
  unhold = true;
 otherwise,
  unhold = false;
end;

xind = [1; xind; length(x)];

%run through each of the categories for x values
for i = 2:length(xind),
    y1 = y(xind(i-1):xind(i)-1);
    gp1 = gp(xind(i-1):xind(i)-1);
    x1 = x(xind(i-1));
    for j = 1:length(gps),
        ox = opt.xoffset(j);
        drawbox(x1+ox,y1(gp1 == gps(j)),opt.color(j,:),...
                opt.width,opt.fill,opt.notch,opt.outliers,opt.marker);
        hold on;
    end;
end;

if (unhold),
    hold off;
end;



% *************************************************
function drawbox(x,y,col,w,fill,notch,outliers,mark)
whis = 1.5;
if (isempty(y)),
    plot(x(1),NaN);
    return;
end;
med = median(y);
q1 = prctile(y,25);
q3 = prctile(y,75);
% find the extreme values (to determine where whiskers appear)
if (outliers),
    vhi = q3+whis*(q3-q1);
    upadj = max(y(y<=vhi));
    if (isempty(upadj)), 
        upadj = q3; 
    end;
    vlo = q1-whis*(q3-q1);
    loadj = min(y(y>=vlo));
    if (isempty(loadj)), 
        loadj = q1; 
    end;
    extreme = y((y < loadj) | (y > upadj));
else
    upadj = max(y);
    loadj = min(y);
    extreme = [];
end;
n1 = med + 1.57*(q3-q1)/sqrt(length(y));
n2 = med - 1.57*(q3-q1)/sqrt(length(y));
if n1>q3, 
    n1 = q3; 
end;
if n2<q1, 
    n2 = q1; 
end;
if (notch),
    boxx = [-0.25 -0.5 -0.5 0.5 0.5 0.25 0.5 0.5 -0.5 -0.5 -0.25]*w + x;
    boxy = [med n1 q3 q3 n1 med n2 q1 q1 n2 med];
    medx = [-0.25 0.25]*w + x;
    medy = [med med];
else
    boxx = [-0.5 -0.5 0.5 0.5 -0.5]*w + x;
    boxy = [q1 q3 q3 q1 q1];
    medx = [-0.5 0.5]*w + x;
    medy = [med med];
end;
whx = [x x NaN x x];
why = [upadj q3 NaN q1 loadj];
if (~fill),
    plot(boxx,boxy,'-',medx,medy,'-',whx,why,'-','Color',col);
    if (~isempty(extreme)),
        addplot(x,extreme,mark,'Color',col);
    end;
else
    patch(boxx,boxy,col,'EdgeColor','none');
    addplot(medx,medy,'Color','w');
    addplot(whx,why,'-','Color',col);
    if (~isempty(extreme)),
        addplot(x,extreme,mark,'Color',col);
    end;
end;
