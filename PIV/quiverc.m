function varargout = quiverc(varargin)
% function quiverc(x,y,u,v,...)
%	or     quiverc(u,v,...)
%	or	   quiverc(...,col,...)
%
% Plots a vector field at coordinates (x,y) with vector length (u,v).  Or
% plots on an evenly spaced grid the vector matrices (u,v).  Can color the
% vectors according to another matrix, col, which must be the same size as u
% and v.  You can also pass a single character color as in other plotting
% routines.  Currently does not handle RGB colors correctly.  Automatically
% scale the vectors so that the 95% percentile vector length corresponds to
% the average spacing between vectors.  Vector head size scales with length,
% also.
%
% Takes a lot of options to customize the display of the vectors.  Options
% can come in any order, and usually take a numeric argument after the
% option name.  All options have a long name and a short name shown in
% parenthesis.
%
% 'ScaleFactor' ('s'),sf - Multiplies the autoscale length by sf.  
%        Default: 1
% 'RelScale' ('rs'),rs - Changes which percentile corresponds to the 
%        maximum(ish) length.  Percentiles are from 0 to 1.  Default: 0.95
% 'AbsScale' ('as'),scale - Absolute scaling value.  Vector display lengths 
%        are created by multiplying the vector length by this value.
% 'ScaleRange' ('sr'),[lo hi] - Changes the range over which vector display
%        length scales by true length.  Expressed in percentiles from 0 to 
%        1.  Lengths below lo are shown as lo in length, and lengths above 
%        hi are only hi long.
% 'AbsScaleRange' ('asr'),[lo hi] - Same as ScaleRange, but expressed in 
%        absolute length, not percentiles.
% 'HeadSize' ('hs'),hsz - Changes the head size as a percentage from 0 to 
%        1 of the total length of the vector.  Default: 0.4.
% 'HeadRange' ('hr'),[hlo hhi] - Like 'ScaleRange', but only affects the 
%        scaling of the head.  A typical use is to have small heads remain 
%        even when the vectors get really small, so that you can still see
%        the direction: 'HeadRange',[0.2 1], which means that all vectors 
%        shorter than the 20th percentile will still have a head.
% 'AbsHeadRange' ('ahr'),[hlo hhi] - Same as 'HeadRange', but takes absolute
%        lengths, not percentiles.
% 'Truncate' ('t'),maxlen - Old option provided for backwards compatibility.
%        Cuts off scaling at maxlen percentile.  Use ScaleRange instead.
% 'Show' ('sh'),fac - Shows fac times fewer vectors (fac < 1).  Mostly use
%        multiples of 1/4, which corresponds to half the number of vectors in
%        both x and y.
% 'NoHeads' ('nh') - Doesn't show vector heads.
% 'NoTails' ('nt') - Doesn't show vector tails.
% 'CorrectAspectRatio' ('car') - Draw the arrow heads so that they look
%        right, even when the axes have different scales
% 'Polar' ('pol') - Vectors are specified in length (u) and angle (v, in
%        radians)
%
% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2010, Eric Tytell

if (isstruct(varargin{1})),
    opts = varargin{1};
    x = opts.x;
    y = opts.y;
    z = opts.z;
    u = opts.u;
    v = opts.v;
    w = opts.w;
    col = opts.Color;
else
    [x,y,z,u,v,w,col,opts] = parseQuivercParameters(varargin);
end;

if (~isfield(opts,'Axes') || isempty(opts.Axes) || (opts.Axes == -1)),
    opts.Axes = gca;
end;

if (strcmp(get(gcf,'Renderer'),'painters')),
  isPainters = true;
else
  isPainters = false;
end;

if (~isPainters),
  if (ischar(col)),
    switch lower(col),
     case {'y','yellow'},
      col = [1 1 0];
     case {'m','magenta'},
      col = [1 0 1];
     case {'c','cyan'},
      col = [0 1 1];
     case {'r','red'},
      col = [1 0 0];
     case {'g','green'},
      col = [0 1 0];
     case {'b','blue'},
      col = [0 0 1];
     case {'w','white'},
      col = [1 1 1];
     case {'k','black'},
      col = [0 0 0];
    end;
  end;
end;

switch get(opts.Axes,'NextPlot'),
 case 'replace',
  cla(opts.Axes,'reset');
 case 'replacechildren',
  cla(opts.Axes);
end;

% definition of the arrow shape
if (~opts.is3D),
    headx = [0.85 	-0.15	0    -0.15    0.85]';
    heady = [0	-0.35	0    0.35     0]';

    if (opts.Angle),
        tailx = [-0.5 0.5]';
        taily = [0 0]';
    else
        tailx = [0		1]';
        taily = [0		0]';
    end;
else
    tailx = [0  1]';
    taily = [0  0]';

    if (opts.HeadFacets3d == 1),
        [~,plane] = min([range(x(:)) range(y(:)) range(z(:))]);
        switch plane,
         case 3,                        % xy plane
          headx = [0.85 -0.15	0    -0.15    0.85]';
          heady = [0	-0.35	0    0.35     0]';
          headz = [0      0      0     0    0]';
         case 2,                        % xz plane
          headx = [0.85 -0.15	0    -0.15    0.85]';
          heady = [0      0      0     0    0]';
          headz = [0	-0.35	0    0.35     0]';
         case 1,                        % yz plane
          headx = [0      0      0     0    0]';
          heady = [0.85 -0.15	0    -0.15    0.85]';
          headz = [0	-0.35	0    0.35     0]';
        end;
    elseif (opts.HeadFacets3d == 2),
        headx = repmat([0.85 -0.15 0]',[1 4]);
        heady = [0 -0.35 0; 0 0.35 0; 0 0 0; 0 0 0]';
        headz = [0 0 0; 0 0 0; 0 -0.35 0; 0 0.35 0]';
    else
        headx1 = repmat([0.85 0 0]',[1 opts.HeadFacets3d]);
        hy = repmat([0 -0.35 -0.35]',[1 opts.HeadFacets3d]);
        s = 0.35*tan(pi/opts.HeadFacets3d);
        %order of z coordinates is arranged so that the normal vector
        %points out
        hz = repmat([0 s -s]',[1 opts.HeadFacets3d]);
        theta = repmat(linspace(0,2*pi,opts.HeadFacets3d+1),[3 1]);
        theta = theta(:,1:end-1);

        heady1 = hy.*cos(theta) - hz.*sin(theta);
        headz1 = hy.*sin(theta) + hz.*cos(theta);
    
        headx2 = repmat([0 0 0.2]',[1 opts.HeadFacets3d]);
        hy = repmat([-0.35 -0.35 0]',[1 opts.HeadFacets3d]);
        s = 0.35*tan(pi/opts.HeadFacets3d);
        %different order here, but also so that the normal points out
        hz = repmat([-s s 0]',[1 opts.HeadFacets3d]);
        theta = repmat(linspace(0,2*pi,opts.HeadFacets3d+1),[3 1]);
        theta = theta(:,1:end-1);

        heady2 = hy.*cos(theta) - hz.*sin(theta);
        headz2 = hy.*sin(theta) + hz.*cos(theta);
        
        headx = [headx1 headx2];
        heady = [heady1 heady2];
        headz = [headz1 headz2];
    end;
end;

tmult = size(tailx,2);
hmult = size(headx,2);

if (~opts.is3D),
    z = zeros(size(x));
    w = zeros(size(u));
end;

% first get the spacing between points in x and y
dx = range(x(:));
dy = range(y(:));
dz = range(z(:));

%get aspect ratio correction factors
xlm = get(opts.Axes,'XLim');
ylm = get(opts.Axes,'YLim');
if (opts.is3D),
    zlm = get(opts.Axes,'ZLim');
    AR = [diff(xlm) diff(ylm) diff(zlm)];
else
    AR = [diff(xlm) diff(ylm)];
end;
AR = AR / mean(AR);

% reduce the number of vectors according to the Show option
if (numel(x) == length(x)),
% random vector positions
  N = length(x);
  
  k = round(linspace(1,N,opts.Show*N));
  x = x(k);
  y = y(k);
  z = z(k);
  u = u(k);
  v = v(k);
  w = w(k);

  if (opts.isColMatrix),
    col = col(k);
  end;

  if (~opts.is3D),
      m = sqrt(length(x));
      n = m;
      p = 0;
  else
      m = length(x)^(1/3);
      n = m;
      p = m;
  end;
else
    % gridded vector positions
    [m,n,p] = size(x);
    
    if (length(opts.Show) == 1),
        if (opts.is3D),
            shw = opts.Show^(1/3);
        else
            shw = sqrt(opts.Show);
        end;
        shw = shw .* [1 1 1];
    elseif (length(opts.Show) == 2),
        shw = [makerow(opts.Show) 1];
    else
        shw = opts.Show;
    end;

    if (mod(1/shw(1),1) == 0),
        step = floor(1/shw(1));
        i = 1:step:m;
        i = i + floor((m-i(end))/2);
    else
        i = round(linspace(1,m,shw(1)*m));
    end;
    if (mod(1/shw(2),1) == 0),
        step = floor(1/shw(2));
        j = 1:step:n;
        j = j + floor((n-j(end))/2);
    else
        j = round(linspace(1,n,shw(2)*n));
    end;
    if (mod(1/shw(3),1) == 0),
        step = floor(1/shw(3));
        k = 1:step:p;
        k = k + floor((p-k(end))/2);
    else
        k = round(linspace(1,p,ceil(shw(3)*p)));
    end;

    x = x(i,j,k);
    y = y(i,j,k);
    z = z(i,j,k);
    u = u(i,j,k);
    v = v(i,j,k);
    w = w(i,j,k);

    if (opts.isColMatrix),
        col = col(i,j,k);
    end;
    
    [m,n,p] = size(x);
    if (~opts.is3D),
        p = 0;
    end;
end;

%%%%%% Set up scaling

% weight spacing by the number of points in each direction
d = (dx + dy + dz)/(n+m+p);
if ((dx == 0) && (dy == 0) && (dz == 0) && isempty(opts.AbsScale)),
  opts.AbsScale = 1;
end;

% Remove any NaNs or zeros
good = isfinite(x) & isfinite(y) & isfinite(z) & ...
         isfinite(u) & isfinite(v) & isfinite(w) & ...
         ((u ~= 0) | (v ~= 0) | (w ~= 0));
if (all(~good(:))),
  warning('quiverc:NoVectors','No vectors to plot.');
  if (nargout == 1),
    varargout = {-1};
  end;
  return;
end;
x = makerow(x(good));
y = makerow(y(good));
z = makerow(z(good));
u = makerow(u(good));
v = makerow(v(good));
w = makerow(w(good));
if (opts.isColMatrix),
  col = makerow(col(good));
end;
N = sum(good(:));

% Set up the vectors so that they have a length (len) and a direction (u,v)
if (opts.Polar),
    len = u;
    ang = v;
    if (any(w ~= 0)),
        error('Cannot plot polar vectors in 3D');
    end;
    u = cos(ang);
    v = sin(ang);
else
    len = sqrt(u.^2 + v.^2 + w.^2);
    u = u./len;
    v = v./len;
    w = w./len;
    if (opts.Angle),
      len = ones(size(len));
    end;
end;

slen = len;

% Define scale range
if (~isempty(opts.AbsScaleRange)),
  isShort = len < opts.AbsScaleRange(1);
  isLong = len > opts.AbsScaleRange(2);
  slen(isShort) = min(len(len >= opts.AbsScaleRange(1)));
  slen(isLong) = max(len(len <= opts.AbsScaleRange(2)));

  len1 = opts.AbsScaleRange(1);
  len2 = opts.AbsScaleRange(2);
else
  short = prctile(len(:),opts.ScaleRange(1)*100);
  long = prctile(len(:),opts.ScaleRange(2)*100);
  isShort = len < short;
  isLong = len > long;
  slen(isShort) = short;
  slen(isLong) = long;

  len1 = short;
  len2 = long;
end;

% Autoscale
if (~isempty(opts.AbsScale)),
  opts.Scale = opts.AbsScale * opts.ScaleFactor;
elseif (~opts.ScaleAsPrevious || ~isfield(opts,'Scale')),
  islen = (len >= len1) & (len <= len2);
  long = prctile(len(islen),opts.RelScale*100);
  opts.Scale = d/long * opts.ScaleFactor;
end;

headlen = slen*opts.HeadSize;
if (~isempty(opts.AbsHeadRange)),
  isShortHead = len < opts.AbsHeadRange(1);
  isLongHead = len > opts.AbsHeadRange(2);
  headlen(isShortHead) = opts.AbsHeadRange(1)*opts.HeadSize;
  headlen(isLongHead) = opts.AbsHeadRange(2)*opts.HeadSize;
else
  short = prctile(len(:),opts.HeadRange(1)*100);
  long = prctile(len(:),opts.HeadRange(2)*100);
  
  isShortHead = len < short;
  isLongHead = len > long;
  headlen(isShortHead) = short*opts.HeadSize;
  headlen(isLongHead) = long*opts.HeadSize;

  if (opts.HeadRange > 0),
    opts.AbsHeadRange(1) = short;
  else
    opts.AbsHeadRange(1) = 0;
  end;
  if (opts.HeadRange < 1),
    opts.AbsHeadRange(2) = long;
  else
    opts.AbsHeadRange(2) = Inf;
  end;
end;

taillen = slen - headlen;
taillen(taillen < 0) = 0;

if (opts.Scale ~= 0),
  taillen = taillen*opts.Scale;
  headlen = headlen*opts.Scale;
else
  taillen = repmat(d*(1-opts.HeadSize),size(len));
  headlen = repmat(d*opts.HeadSize,size(len));
end;

%%%% Construct vectors

if (opts.NoTails),
    headlen = headlen+taillen;
    taillen = zeros(size(taillen));
elseif (opts.NoHeads),
    taillen = headlen+taillen;
    headlen = zeros(size(headlen));
end;

if (~opts.NoTails),
  tx = repmat(tailx,[1 N]);
  ty = repmat(taily,[1 N]);
  
  len = flatten(repmat(taillen, [tmult 1]));
  nonZeroTail = len > 0;

  tx = tx.*repmat(len',[size(tx,1) 1]);
  ty = ty.*repmat(len',[size(tx,1) 1]);
  if (opts.is3D),
      tz = zeros(size(tx));
  end;
  nt = size(tx,1);
end;

if (~opts.NoHeads),
  hx = repmat(headx,[1 N]);
  hy = repmat(heady,[1 N]);

  len = flatten(repmat(headlen,[hmult 1]));
  hx = hx.*repmat(len',[size(hx,1) 1]);
  hy = hy.*repmat(len',[size(hx,1) 1]);
  if (opts.is3D),
      hz = repmat(headz,[1 N]);
      hz = hz.*repmat(len',[size(hx,1) 1]);
  end;
  nh = size(hx,1);
end;

if (~opts.NoTails),
    uu = repmat(flatten(repmat(u,[tmult 1]))',[nt 1]);
    vv = repmat(flatten(repmat(v,[tmult 1]))',[nt 1]);
    xx = repmat(flatten(repmat(x,[tmult 1]))',[nt 1]);
    yy = repmat(flatten(repmat(y,[tmult 1]))',[nt 1]);
    
    if (opts.is3D),
        zz = repmat(flatten(repmat(z,[tmult 1]))',[nt 1]);
        ww = repmat(flatten(repmat(w,[tmult 1]))',[nt 1]);
        aa = sqrt(uu.^2 + vv.^2);
        
        tailptx = uu.*tx - vv./aa.*ty - ww.*uu./aa.*tz;
        tailpty = vv.*tx + uu./aa.*ty - ww.*vv./aa.*tz;
        tailptz = ww.*tx + aa.*tz;
    else
        tailptx = uu.*tx - vv.*ty;
        tailpty = vv.*tx + uu.*ty;
    end;

    if (opts.Polar && opts.CorrectAspectRatio),   
        tailptx = tailptx * AR(1);
        tailpty = tailpty * AR(2);
        if (opts.is3D),
            tailptz = tailptz * AR(3);
        end;
    end;
    tailptx = tailptx + xx;
    tailpty = tailpty + yy;
    if (opts.is3D),
        tailptz = tailptz + zz;
    end;
end;
    
if (~opts.NoHeads),
    uu = repmat(flatten(repmat(u,[hmult 1]))',[nh 1]);
    vv = repmat(flatten(repmat(v,[hmult 1]))',[nh 1]);
    
    if (~opts.NoTails)
        %add tail x and y to xx and yy
        tx1 = tailptx(nt,1:tmult:end);
        ty1 = tailpty(nt,1:tmult:end);
        
        xx = repmat(flatten(repmat(tx1,[hmult 1]))',[nh 1]);
        yy = repmat(flatten(repmat(ty1,[hmult 1]))',[nh 1]);
    else
        xx = repmat(flatten(repmat(x,[hmult 1]))',[nh 1]);
        yy = repmat(flatten(repmat(y,[hmult 1]))',[nh 1]);
    end;
    
    if (opts.is3D),
        zz = repmat(flatten(repmat(z,[hmult 1]))',[nh 1]);
        ww = repmat(flatten(repmat(w,[hmult 1]))',[nh 1]);
        aa = sqrt(uu.^2 + vv.^2);
        
        headptx = uu.*hx - vv./aa.*hy - ww.*uu./aa.*hz;
        headpty = vv.*hx + uu./aa.*hy - ww.*vv./aa.*hz;
        headptz = ww.*hx + aa.*hz;
    else
        headptx = uu.*hx - vv.*hy;
        headpty = vv.*hx + uu.*hy;
    end;
    
    if (opts.CorrectAspectRatio),
        headptx = headptx * AR(1);
        headpty = headpty * AR(2);
        if (opts.is3D),
            headptz = headptz * AR(3);
        end;
    end;

    headptx = headptx + xx;
    headpty = headpty + yy;
    if (opts.is3D),
        headptz = headptz + zz;
    end;
end;

if (opts.isColMatrix),
    if (~opts.NoHeads),
        headcol = repmat(flatten(repmat(col,[hmult 1]))',[nh 1]);
    end;
    if (~opts.NoTails),
        tailcol = repmat(flatten(repmat(col,[tmult 1]))',[nt 1]);
    end;
end;

hind = 1;
if (~opts.NoHeads),
    if (~opts.isColMatrix),
        headcol = col;
    end;

    if (opts.is3D),
        h(hind) = patch(headptx,headpty,headptz, headcol,...
                        'EdgeColor','none', 'Parent',opts.Axes);
    else
        h(hind) = patch(headptx,headpty, headcol,...
                          'EdgeColor','none', 'Parent',opts.Axes);
    end;
    hind = hind+1;
end;
if (~opts.NoTails),
  if (opts.isColMatrix),
      if (opts.is3D),
          h(hind) = ...
              patch(tailptx(:,nonZeroTail),tailpty(:,nonZeroTail),...
                    tailptz(:,nonZeroTail),'b',...
                    'FaceColor','none',...
                    'FaceVertexCData', flatten(tailcol(:,nonZeroTail)),...
                    'EdgeColor','flat','LineWidth',opts.LineWidth, ...
                    'Parent',opts.Axes);
      else
          h(hind) = ...
              patch(tailptx(:,nonZeroTail),tailpty(:,nonZeroTail),'b',...
                    'FaceColor','none',...
                    'FaceVertexCData', flatten(tailcol(:,nonZeroTail)),...
                    'EdgeColor','flat','LineWidth',opts.LineWidth,...
                    'Parent',opts.Axes);
      end;
  else
    if (opts.NoHeads),
      xx = flatten([tailptx(:,nonZeroTail); ...
                        repmat(NaN,[1 sum(nonZeroTail)])]);
      yy = flatten([tailpty(:,nonZeroTail); ...
                        repmat(NaN,[1 sum(nonZeroTail)])]);
      if (opts.is3D),
          zz = flatten([tailptz(:,nonZeroTail); ...
                        repmat(NaN,[1 sum(nonZeroTail)])]);
          hold on;
          h(hind) = plot3(opts.Axes, xx,yy,zz, 'Color',col,'LineWidth',opts.LineWidth);
      else
          h(hind) = addplot(opts.Axes, xx,yy, 'Color',col,'LineWidth',opts.LineWidth);
      end;
    else
        if (any(nonZeroTail)),
            if (opts.is3D),
                h(hind) = ...
                    patch(tailptx(:,nonZeroTail),tailpty(:,nonZeroTail),...
                    tailptz(:,nonZeroTail),...
                    'k','FaceColor','none','EdgeColor',col,...
                    'LineWidth',opts.LineWidth, 'Parent',opts.Axes);
            else
                h(hind) = ...
                    patch(tailptx(:,nonZeroTail),tailpty(:,nonZeroTail),...
                    'k','FaceColor','none','EdgeColor',col,...
                    'LineWidth',opts.LineWidth, 'Parent',opts.Axes);
            end;
        end;
    end;
  end;
end;

set(h,'Tag','quiverc', 'UserData',opts);

if (opts.is3D && strcmp(get(opts.Axes,'NextPlot'),'replace')),
    view([-37.5 30]);
    camlight;
end;

if (nargout > 0),
  varargout{1} = h;
  if (nargout == 2),
      varargout{2} = opts;
  end;
end;

