function pivplot(x,y,u,v,varargin)

% set up the data matrix
D.x = x;
D.y = y;
D.u = u;
D.v = v;

% defaults
D.imbg = [];
D.bkgnd = [];
D.color = [];
ctr = 'auto';
rng = [];
step = [];
axistype = 'xy';
D.iscolorbar = 1;

D.quivercopts = {};

% parse optional parameters
opts = varargin;
i = 1;
while (i <= length(opts)),
    switch lower(opts{i}),
     case {'image','im'},
      D.imbg = opts{i+1};
      i = i+2;

     case {'background','bkgnd','bg'},
      D.bkgnd = opts{i+1};
      i = i+2;
     case {'color','col'},
      D.color = opts{i+1};
      i = i+2;

     case {'center','ctr'},
      if (isnumeric(opts{i+1}) | ...
          (ischar(opts{i+1}) & ...
           ~isempty(strmatch(opts{i+1}, {'auto','none','median'},...
                             'exact')))),
          ctr = opts{i+1};
          i = i+2;
      else
          ctr = 'auto';
          i = i+1;
      end;
     case {'range','rng'},
      rng = opts{i+1};
      i = i+2;

     case 'colormap',
      cmap = opts{i+1};
      i = i+2;
     case 'cmapstep',
      cmapstep = opts{i+1};
      i = i+2;

     case {'xy','ij'},
      axistype = opts{i};
      i = i+1;

    % quiverc options
     case {'show','sh','truncate','t','linewidth','lw','headsize','hs',...
           'headrange','hr','absheadrange','ahr','scalefactor','s',...
           'absscale','as'},
      D.quivercopts(end+1:end+2) = opts(i:i+1);
      i = i+2;

     otherwise,
      if (ischar(opts{i})),
          error('Unrecognized option %s.');
      else
          error('Unrecognized option number %d.',i+4);
      end;
    end;
end;

% check sizes on everything
if (any(size(D.u) ~= size(D.v))),
    error('u must be the same size as v');
end;
if (any(size(D.x) == 1) & (numel(D.x) ~= size(D.u,2))),
    error(['If x is a vector, it must have the same number of colums as ' ...
           'the vector matrices.']);
elseif (any([size(D.x,1) size(D.x,2)] ~= [size(D.u,1) size(D.u,2)]) | ...
        ((size(D.x,3) ~= 1) & (size(D.x,3) ~= size(D.u,3)))),
    error('x must be the same size as the vector matrices.');
end;
if (any(size(D.y) == 1) & (numel(D.y) ~= size(D.u,1))),
    error(['If y is a vector, it must have the same number of colums as ' ...
           'the vector matrices.']);
elseif (any([size(D.y,1) size(D.y,2)] ~= [size(D.u,1) size(D.u,2)]) | ...
        ((size(D.y,3) ~= 1) & (size(D.y,3) ~= size(D.u,3)))),
    error('y must be the same size as the vector matrices.');
end;

if (ischar(D.bkgnd)),
    D.bkgnd = generateScalar(D.bkgnd, D.x,D.y,D.u,D.v);
end;
if (ischar(D.color)),
    D.color = generateScalar(D.color, D.x,D.y,D.u,D.v);
end;

if (~isempty(D.bkgnd) & ...
    ((size(D.bkgnd,1) ~= size(D.u,1)) | ...
     (size(D.bkgnd,2) ~= size(D.u,2)) | ...
     ((size(D.bkgnd,3) ~= 1) & (size(D.bkgnd,3) ~= size(D.u,3))))),
    error(['Background matrix must be the same size as the vector ' ...
           'matrix.']);
end;
if (~isempty(D.color) & ...
    ((size(D.color,1) ~= size(D.u,1)) | ...
     (size(D.color,2) ~= size(D.u,2)) | ...
     ((size(D.color,3) ~= 1) & (size(D.color,3) ~= size(D.u,3))))),
    error(['Color matrix must be the same size as the vector ' ...
           'matrix.']);
end;

if (~isempty(D.bkgnd) & ~isempty(D.color)),
    error('Cannot display both a background and colored vectors (yet)');
end;

switch get(gca,'NextPlot'),
 case 'replace',
  cla reset;
 case 'replacechildren',
  cla;
end;
        
% prepare range
if (~isempty(rng) & (~isempty(D.bkgnd) | ~isempty(D.color))),
    if (numel(rng) == 2),
        D.caxlim = rng;
    else
        if (~isempty(D.bkgnd)),
            it = D.bkgnd;
        else
            it = D.color;
        end;

        if (strcmp(ctr,'auto')),
            ctr = nanmean(it(:));
        elseif (strcmp(ctr,'median')),
            ctr = nanmedian(it(:));
        end;
        
        if (~isempty(rng)),
            stdev = nanstd(it(:));
            D.caxlim = ctr + rng*[-stdev stdev];
        else
            m = max(abs(it(:)-ctr));
            D.caxlim = ctr + [-m m];
        end;

        if (numel(cmapstep) == 1),
            symmetricColormap(cmap,cmapstep/(2*rng));
        else
            smallstep = min(cmapstep);
            cmapctr = cmapstep(1)/smallstep;
            cmapstep = (rng-cmapstep(1)/2)/smallstep - 1;

            symmetricColormap(cmap,cmapctr,cmapstep);
        end;
    end;
end;

D.hBkgnd = -1;
D.hImage = -1;
D.hVectors = -1;
D.fr = 1;
D.nfr = size(D.u,3);

if (size(D.u,3) > 1),
    D.ax = gca;
    D.fig = gcf;
    set(D.fig,'KeyPressFcn',@pivplotKeyPress, 'UserData','');
    D.axlim = [];

    guidata(D.fig,D);

    D = dopivplot(D,1);

    D.axlim(1:2) = get(D.ax,'XLim');
    D.axlim(3:4) = get(D.ax,'YLim');
    axis('equal',axistype);
    axis(D.axlim);

    guidata(D.fig,D);

    uiwait(D.fig);
    set(D.fig,'KeyPressFcn','');
else
    dopivplot(D,1);
end;



% ------------------------------------------------------------------
function D = dopivplot(D,fr)

if (~isempty(D.bkgnd)),
    if (size(D.bkgnd,3) == 1),
        bgfr = 1;
    else
        bgfr = fr;
    end;

    if (ishandle(D.hBkgnd)),
        set(D.hBkgnd,'CData',D.bkgnd(:,:,bgfr));
    else
        D.hBkgnd = imagesc(D.x(1,:),D.y(:,1),D.bkgnd(:,:,bgfr));
    end;
end;

if (~isempty(D.color)),
    if (size(D.color,3) == 1),
        qcol = {D.color};
    else
        qcol = {D.color(:,:,fr)};
    end;
else
    qcol = {'k'};
end;

if (ishandle(D.hVectors)),
    delete(D.hVectors);
end;

hold on;
D.hVectors = ...
    quiverc(D.x,D.y,D.u(:,:,fr),D.v(:,:,fr),qcol{:},D.quivercopts{:});
hold off;
caxis(D.caxlim);

if (D.iscolorbar),
    colorbar;
end;

% ------------------------------------------------------------------
function s = generateScalar(what, x,y,u,v)

switch lower(what),
 case 'vorticity',
  s = vorticity(x,y,u,v);

 otherwise,
  error('Unknown scalar type %s.',what);
end;

% ------------------------------------------------------------------
function pivplotKeyPress(obj,eventdata)

D = guidata(obj);
c = get(D.fig, 'CurrentCharacter');

draw = 0;
switch lower(c),
 case char(28),                         % left arrow
  if (D.fr > 1),
    D.fr = D.fr-1;
  end;
 case {char(29),char(13),' '},          % right arrow, return
  if (D.fr < D.nfr),
    D.fr = D.fr+1;
  end;
 case 'g',                              % go to frame
  fr = inputdlg('Frame?','Go to frame');
  fr = str2num(fr{1});
  if ((fr >= 1) & (fr <= D.nfr)),
    D.fr = fr;
  end;
 case 'q',                              % quit
  uiresume;
end;

D = dopivplot(D,D.fr);
guidata(obj,D);





            
        












