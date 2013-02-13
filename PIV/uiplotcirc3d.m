function uiplotCirc3D(vxr, varargin)
% function uiplotCirc3D(vxr,[vxrGood], x,y,t,u,v,[dcev], ...)
%          uiplotCirc3D(..., 'flowvel',flowvel)
%             Defines the flow velocity to convert t into an artificial
%             z position.  Without it, just uses t directly.
%          uiplotCirc3D(..., 'body',bodyx,bodyy)
%             Matrix of body position.  Can have any number of rows, but
%             must have the same number of columns as u has frames, or else
%             you must pass in pivframes 
%          uiplotCirc3D(..., 'pivframes',pivframes)
%             Defines the correspondence between t (and bodyx and bodyy,
%             if they're passed in) and frames in u and v.
%          uiplotCirc3D(..., 'frames',showframes)
%             Range in frames to show.
%          uiplotCirc3D(..., 'range',showrange)
%             Range in t to show.
%          uiplotCirc3D(..., 'showvec',showvec)
%             Sets the percentage of vectors to show.  Usually 0.25.
%          uiplotCirc3D(..., 'avi',imx,imy,aviname,aviframes)
%             Shows an avi in the background.  imx are the coordinates of
%             the image edges; aviname is the name of the file; and 
%             aviframes gives the frames in the avi that correspond to 
%             the pages of u and v.
%
% Click on a vortex in either window to select it and display it 
% highlighted in the other window.  Hit delete to delete the vortex 
% (although it remains in the vxr matrix).
% You can double click on vortices in the 2D window to change their 
% names.  Also right clicking gives a context menu that lets you change
% the start and end frame of vortices (which changes the vxr structure
% irreversibly at the moment).
%
% Saves the vxr structure and a logical array identifying which vectors 
% were deleted back into the base namespace.
%
% Current issues: 
% * For some reason, maximizing either window tends to make Matlab's OpenGL 
%   renderer unhappy.  If you just resize them manually to take up the whole 
%   screen, it seems OK, though.
% * In principle, you can double click in either window to change the vortex
%   name.  But if you've got the camera toolbar enabled in the 3D window,
%   double clicking severely confuses Matlab and ends up permanently
%   disabling the camera toolbar.  Soln: just double click in the 2D window.

%save the name of the vortex matrix
varnames{1} = inputname(1);

% check whether they gave us the vxrGood array
pnum = 1;
if ((ndims(vxr) == ndims(varargin{pnum})) & ...
    (length(varargin{pnum}) == length(vxr)) & ...
    islogical(varargin{pnum})),
    vxrGood = varargin{pnum};
    pnum = pnum+1;
    varnames{2} = inputname(2);
else
    vxrGood = logical(ones(size(vxr)));
    varnames{2} = '';
end;

% now get x,y,z,u, and v
if (length(varargin) < pnum+4),
    error('Not enough arguments to uiplotCirc3D');
end;
x = varargin{pnum};
y = varargin{pnum+1};
z = varargin{pnum+2};
z0 = z;
u = varargin{pnum+3};
v = varargin{pnum+4};
pnum = pnum+5;

%check whether they passed a dcev
if ((length(varargin) >= pnum) & (ndims(varargin{pnum}) == ndims(u)) & ...
    all(size(varargin{pnum}) == size(u)) & isnumeric(varargin{pnum})),
    dcev = varargin{pnum};
    pnum = pnum+1;
end;

%assume there isn't a body outline, unless we get one below
isBody = false;
isAvi = false;
showVec = 0.25;                         % show 25% of the vectors

% check on really optional things
opts = varargin(pnum:end);
i = 1;
while (i <= length(opts)),
    %error if we don't have a string identifying the option
    if (~ischar(opts{i})),
        error('Unrecognized option %d.',i+pnum+6);
    end;

    switch lower(opts{i}),
     case {'flow','flowvel'},
      %the z var is actually a time - convert to a position by using the
      %flow velocity
      flowvel = opts{i+1};
      z = z0*flowvel;
      i = i+2;

     case 'body',
      %they passed in an outline of the body for each frame
      bodyx = opts{i+1};
      bodyy = opts{i+2};
      isBody = true;
      i = i+3;

     case {'showframes','frames'},
      %frames to show
      showFrames = opts{i+1};
      i = i+2;

     case 'pivframes',
      %frames that correspond between t and u
      pivframes = opts{i+1};
      i = i+2;

     case {'showrange','range'},
      %range in z0 to show
      showRange = opts{i+1};
      i = i+2;

     case 'showvec',
      %percentage of the vectors to show in the 2D plot
      showVec = opts{i+1};
      i = i+2;

     case 'avi',
      %show avi frames in the background
      %avix and aviy are the two extreme positions of the edges of the
      %avi.  If they're more than 2 element vectors, we're just going to
      %take the ends, anyway
      avix = opts{i+1};
      aviy = opts{i+2};
      aviname = opts{i+3};
      aviframes = opts{i+4};
      isAvi = true;
      i = i+5;

    end;
end;
      
%check sizes
if (length(vxr) ~= length(vxrGood)),
    error('vxrGood must be a logical array the same size as vxr');
end;
if (any(size(x) ~= size(y))),
    error('x and y must be the same sizes');
end;
if (any([size(x,1) size(x,2)] ~= [size(u,1) size(u,2)])),
    error(['The first two dimensions of x and y must be the same as u and ' ...
           'v.']);
end;
if (any(size(u) ~= size(v))),
    error('u and v must be the same size');
end;
if (exist('dcev') & (any(size(dcev) ~= size(u)))),
    error('dcev must be the same size as u and v');
end;

if (exist('pivframes')),
    npiv = sum(isfinite(pivframes(:)));
    if (size(u,3) ~= npiv),
        error(['pivframes must have the same number of finite elements as '...
               'pages in u, v and dcev']);
    end;
else
    npiv = size(u,3);
end;

if (isBody & (any(size(bodyx) ~= size(bodyy)))),
    error('bodyx and bodyy must be the same size');
end;

if (length(z) ~= npiv),
    if (~exist('pivframes')),
        error('If z has more frames than u, you must pass pivframes');
    end;
    z = z(pivframes(isfinite(pivframes)));
    z0 = z0(pivframes(isfinite(pivframes)));
end;
if (isBody & (size(bodyx,2) ~= npiv)),
    if (~exist('pivframes')),
        error('If bodyx has more frames than u, you must pass pivframes');
    end;
    bodyx = bodyx(:,pivframes(isfinite(pivframes)));
    bodyy = bodyy(:,pivframes(isfinite(pivframes)));
end;
if (isAvi & ~exist(aviname,'file')),
    error('Cannot find avi file %s',aviname');
end;

if (~exist('dcev')),
    %create the dcev if it wasn't passed in
    dcev = disccomplexeig(x,y,u,v);
end;

% restrict to the frame range or z range, if they passed those
if (exist('showRange')),
    if (exist('showFrames')),
        error('Cannot pass both showRange and showFrames');
    end;
    if (numel(showRange) > 2),
        warning(['showRange should only have 2 elements.  Only using ' ...
                 'first and last']);
    end;

    [q,fr1] = min(abs(showRange(1)-z0));
    [q,fr2] = min(abs(showRange(end)-z0));
    showFrames = fr1:fr2;
end;
if (exist('showFrames')),
    if (numel(showFrames) == 2),
        showFrames = showFrames(1):showFrames(end);
    end;

    if (any(showFrames < 1) | any(showFrames > npiv)),
        error('Frames in showFrames are out of range');
    end;

    %grab the appropriate frames
    u = u(:,:,showFrames);
    v = v(:,:,showFrames);
    dcev = dcev(:,:,showFrames);
    z = z(showFrames);
    z0 = z0(showFrames);
    bodyx = bodyx(:,showFrames);
    bodyy = bodyy(:,showFrames);
else
    showFrames = 1:npiv;
end;

fig3D = findobj('Tag','UICirc3DFig');
if (isempty(fig3D)),
    fig3D = figure;
    resetView = true;
else
    resetView = false;
end;
set(fig3D, 'Tag','UICirc3DFig', 'KeyPressFcn',@ucKeyPressFcn, ...
                   'CloseRequestFcn',@ucCloseFig);
fig2D = findobj('Tag','UICirc2DFig');
if (isempty(fig2D)),
    fig2D = figure;
end;
set(fig2D,'Tag','UICirc2DFig', 'KeyPressFcn',@ucKeyPressFcn,...
                   'CloseRequestFcn',@ucCloseFig);

if (min(dcev(:)) < 0),
    dcevRange = [-2000 2000];
    vecCol = 'k';
else
    dcevRange = [0 2000];
    vecCol = 'w';
end;
if (isAvi),
    vecCol = 'y';
end;

%add names to the vxr structure if they aren't already there
if (~isfield(vxr,'name')),
    for i = 1:length(vxr),
        vxr(i).name = '';
    end;
end;

% build the surfaces
figure(fig3D);
if (resetView),
    clf;
else
    cla;
end;
hold on;

hvx3D = -1*ones(size(vxr));
phi = linspace(0,2*pi,17)';

for i = 1:length(vxr),
    %skip vortices that have already been eliminated
    if (~vxrGood(i)),
        continue;
    end;

    inrange = find(ismember(vxr(i).frames,showFrames));
    if (length(inrange) < 2),
        %skip out of range vortices
        continue;
    end;

    nfr = length(inrange);
    %assume for the moment that showFrames is just something like 12:27
    %or whatever
    fr = vxr(i).frames(inrange) - showFrames(1) + 1;

    %get the angle mod pi
    ang1 = mod(vxr(i).angle(inrange),pi);

    %this is the value that phi will have when we're actually at the
    %angle ang1.  It's different because a ~= b
    a = vxr(i).majordiam(inrange)/2;
    b = vxr(i).minordiam(inrange)/2;
    ang2 = atan2(a.*sin(ang1), b.*cos(ang1));

    x1 = zeros(length(phi),nfr);
    y1 = zeros(length(phi),nfr);
    for j = 1:nfr,
        fr1 = inrange(j);

        x1(:,j) = a(j)*cos(phi-ang2(j)).*cos(ang1(j)) - ...
                  b(j)*sin(phi-ang2(j)).*sin(ang1(j)) + ...
                  vxr(i).ctrx(fr1);
        y1(:,j) = a(j)*cos(phi-ang2(j)).*sin(ang1(j)) + ...
                  b(j)*sin(phi-ang2(j)).*cos(ang1(j)) + ...
                  vxr(i).ctry(fr1);
    end;

    z1 = repmat(z(fr),[length(phi) 1]);

    hvx3D(i) = surf(x1,y1,z1,repmat(vxr(i).circmax(inrange),...
                                    [length(phi) 1]),...
                  'EdgeColor','none','UserData',i, ...
                  'ButtonDownFcn',@ucSelectSurface);
    if (isempty(vxr(i).name)),
        alpha(hvx3D(i),0.5);
    end;
end;

if (isBody),
    hbody = surf(bodyx', bodyy', repmat(z,[size(bodyx,1) 1])', ...
                 'FaceColor','k','EdgeColor','k', ...
                 'HitTest','off');
    alpha(hbody,0.5);
end;

if (resetView),
    baseview = [ 0.7814    0.0000    1.7105   -1.2460; ...
                 0.2133    0.9661   -0.4178   -0.3808; ...
                 0.7980   -0.2582   -1.5633    8.5702; ...
                 0         0         0    1.0000];
    view(baseview);
    camup([0 1 0]);                         % set the y axis pointing up
    lighting flat;                          % set the lighting mode
    grid;
    axis equal vis3d tight;
end;
camlight;                               % add a light
axlim = axis;
axlim = axlim(1:4);

% set up guidata structure
data.fig3D = fig3D;
data.fig2D = fig2D;
data.hvx3D = hvx3D;
data.vxr = vxr;
data.vxrGood = vxrGood;
% save the input names of the first two parameters - we'll use them
% as defaults for saving the same parameters back out again
data.varnames = varnames;

data.curFrame = 1;
data.curVx = [];
data.showFrames = showFrames;
data.showVec = showVec;
data.x = x;
data.y = y;
data.z = z;
data.u = u;
data.v = v;
data.vecCol = vecCol;
data.dcev = dcev;
data.dcevRange = dcevRange;
data.axlim = axlim;

if (isBody),
    data.bodyx = bodyx;
    data.bodyy = bodyy;
    data.hbody = [hbody -1];
end;

data.isAvi = isAvi;
if (isAvi),
    data.aviX = avix;
    data.aviY = aviy;
    data.aviName = aviname;
    data.aviFrames = aviframes;
    data.hImage = -1;
end;

data.hDCEV = -1;
data.hVectors = -1;
data.hvx2D = -1*ones(length(vxr),2);
data.hvxlabel = -1*ones(size(vxr));
data.hSlice = -1;

% set up the name edit uicontrol
data.hNameEdit = uicontrol('Parent',fig2D, 'Style','edit','String','',...
                           'Visible','off', 'Position',[20 20 100 20], ...
                           'Callback',@ucNameEditCallback);

%set up the vortex context menu
data.hVortexMenu = uicontextmenu('Visible','off');
uimenu(data.hVortexMenu, 'Label','Set start frame', ...
       'Callback',{@ucMenuCallback,'start'});
uimenu(data.hVortexMenu, 'Label','Set end frame', ...
       'Callback',{@ucMenuCallback,'end'});
uimenu(data.hVortexMenu, 'Label','Split', ...
        'Callback',{@ucMenuCallback,'split'});
% uimenu(data.hVortexMenu, 'Label','Restore', ...
%        'Callback',{@ucMenuCallback,'restore'});
uimenu(data.hVortexMenu, 'Label','Delete', ...
       'Callback',{@ucMenuCallback,'delete'});

guidata(fig3D,data);
guidata(fig2D,fig3D);

figure(fig2D);
cla;
if (dcevRange(1) == 0),
    colormap hot;
else
    colormap default;
end;

data = ucUpdate2D(data);
guidata(fig3D,data);

% -----------------------------------------------------------
function data = ucUpdate2D(data, what)

if ((nargin == 1) | strcmp(what,'all')),
    what = {'background','vectors','body','vortices','slice'};
end;
if (ischar(what)),
    what = {what};
end;

%plot the 2D slice
figure(data.fig2D);
hold on;

fr = data.curFrame;

for i = 1:length(what),
    switch what{i},
     case 'background',
      if (data.isAvi),
          frm = aviread(data.aviName, data.aviFrames(fr));
          %need to have 24bit RGB (3 planes) to get the colormap to work on the
          %DCEV, but not on the image.  If cdata has only one plane, then we
          %need to multiply it to 3 to get it to work
          if (size(frm.cdata,3) ~= 3),
              frm.cdata = repmat(frm.cdata(:,:,1),[1 1 3]);
          end;

          if (ishandle(data.hImage)),
              set(data.hImage, 'CData',frm.cdata);
          else
              data.hImage = image(data.aviX([1 end]),data.aviY([1 end]), ...
                                   frm.cdata);
              axis equal;
          end;
      end;

      if (ishandle(data.hDCEV)),
          set(data.hDCEV, 'CData',data.dcev(:,:,fr));
      else
          data.hDCEV = imagesc(data.x(1,:),data.y(:,1),data.dcev(:,:,fr));
          caxis(data.dcevRange);
          %for some reason, the colorbar can't cope with having the AVI
          %in the background, so only show it if we don't have an image
          if (~data.isAvi),
              colorbar;
          end;
          axis equal;
          axis(data.axlim);
      end;
      
      if (data.isAvi),
          ad = data.dcev(:,:,fr);
          ad(abs(ad) < 100) = 0;
          ad(ad ~= 0) = 0.5;
          
          set(data.hDCEV, 'AlphaData',ad, 'AlphaDataMapping','none');
      end;

     case 'vectors',
      if (ishandle(data.hVectors)),
          delete(data.hVectors);
      end;
      data.hVectors = quiverc(data.x,data.y,data.u(:,:,fr),data.v(:,:,fr),...
                              data.vecCol,'show',data.showVec,'s',0.75);

     case 'body',
      if (isfield(data,'bodyx')),
          if (ishandle(data.hbody(2))),
              set(data.hbody(2), 'XData',data.bodyx(:,fr),...
                                'YData',data.bodyy(:,fr));
          else
              data.hbody(2) = plot(data.bodyx(:,fr),...
                                   data.bodyy(:,fr),'r-');
          end;
      end;

     case 'vortices',
      delete(data.hvx2D(ishandle(data.hvx2D)));
      delete(data.hvxlabel(ishandle(data.hvxlabel)));
      data.hvx2D = -1*ones(length(data.vxr),2);

      figure(data.fig3D);
      hold on;

      goodind = find(data.vxrGood);
      k = findvxrinframe(data.vxr(goodind),data.showFrames(fr));
      k = goodind(k);
      for ii = 1:length(k),
          i = k(ii);

          if (~ishandle(data.hvx3D(i))),
              %only will happen when we've got a vortex that's right at the
              %beginning or end of the frame range that we're showing, such
              %that only 1 or 2 frames of it would be visible.  Then we don't
              %show it in the 3D view because surf gets unhappy trying to
              %render it
              continue;
          end;
          if (data.vxr(i).frames(1) < data.showFrames(1)),
              inrange = find(ismember(data.vxr(i).frames, data.showFrames));
              fr1 = find(data.vxr(i).frames(inrange) == data.showFrames(fr));
          else
              fr1 = find(data.vxr(i).frames == data.showFrames(fr));
          end;

          xd = get(data.hvx3D(i),'XData');
          yd = get(data.hvx3D(i),'YData');
          zd = get(data.hvx3D(i),'ZData');

          figure(data.fig2D);
          data.hvx2D(i,1) = plot(xd(:,fr1),yd(:,fr1),'r-', 'UserData',i, ...
                                 'ButtonDownFcn',@ucSelectSurface);
          if (~isempty(data.vxr(i).name)),
              data.hvxlabel(i) = text(xd(end,fr1),yd(end,fr1), ...
                                      data.vxr(i).name,...
                                      'Color',data.vecCol,...
                                      'ButtonDownFcn',@ucSelectSurface);
          end;

          figure(data.fig3D);
          data.hvx2D(i,2) = plot3(xd(:,fr1),yd(:,fr1),zd(:,fr1),'r-', ...
                                  'HitTest','off');
      end;

      %if the current vortex is visible, highlight it
      if (~isempty(data.curVx) & any(k == data.curVx)),
          set(data.hvx2D(data.curVx), 'Color','y','LineWidth',2);
      end;

     case 'slice',
      figure(data.fig3D);
      [xslice,yslice] = meshgrid(data.axlim(1:2),data.axlim(3:4));
      if (ishandle(data.hSlice)),
          set(data.hSlice,'XData',xslice,'YData',yslice,'ZData',...
                          repmat(data.z(fr),[2 2]));
      else
          data.hSlice = surface(xslice,yslice,...
                                repmat(data.z(fr),[2 2]), ...
                                'FaceColor','r','EdgeColor','r',...
                                'HitTest','off');
          alpha(data.hSlice,0.2);
      end;
    end;
end;

figure(data.fig3D);
hold off;

figure(data.fig2D);
hold off;

% -----------------------------------------------------------
function ucKeyPressFcn(obj,eventdata)

h = guidata(obj);
if (ishandle(h)),
    data = guidata(h);
else
    data = h;
end;

c = get(obj,'CurrentCharacter');

switch lower(c),
 case {char(29),'n'},                   % right arrow or 'n'
  if (data.curFrame < length(data.showFrames)),
      data.curFrame = data.curFrame+1;
  end;
  what = 'all';

 case {char(28),'p'},                    % left arrow or 'p'
  if (data.curFrame > 1),
      data.curFrame = data.curFrame-1;
  end;
  what = 'all';

 case {char(127)},                  % delete
  if (ishandle(data.hvx3D(data.curVx))),
      delete(data.hvx3D(data.curVx));
  end;
  if (ishandle(data.hvx2D(data.curVx))),
      delete(data.hvx2D(data.curVx,:));
  end;
  data.vxrGood(data.curVx) = false;
  data.curVx = [];
  what = 'vortices';

 case 'q',
  ucQuit(data);
  return;
end;

data = ucUpdate2D(data, what);

guidata(data.fig3D,data);

% -----------------------------------------------------------
function ucSelectSurface(obj,eventdata)

h = guidata(obj);
if (ishandle(h)),
    data = guidata(h);
    clicked3D = false;
else
    data = h;
    clicked3D = true;
end;

ind = get(obj,'UserData');

% un-highlight the old surface
set(data.hvx3D(data.curVx), 'EdgeColor','none');
% highlight the 3D surface
set(data.hvx3D(ind), 'EdgeColor','y');

if (clicked3D),
    %find the approximate frame in which we clicked
    c = get(get(obj,'Parent'),'CurrentPoint');
    zclick = sum(c(:,3))/2;
    [q,fr] = min(abs(zclick - data.z));

    %check to make sure we got a frame where our surface is defined
    if (data.showFrames(fr) < data.vxr(ind).frames(1)),
        fr = data.vxr(ind).frames(1) - data.showFrames(1)+1;
    elseif (data.showFrames(fr) > data.vxr(ind).frames(end)),
        fr = data.vxr(ind).frames(end) - data.showFrames(1)+1;
    end;

    data.curVx = ind;
    if (fr ~= data.curFrame),
        data.curFrame = fr;
        data = ucUpdate2D(data, 'all');
    else
        data = ucUpdate2D(data, 'vortices');
    end;
else
    if (all(ishandle(data.hvx2D(data.curVx,:)))),
        set(data.hvx2D(data.curVx,:),'Color','r','LineWidth',0.5);
    end;
    data.curVx = ind;
    set(data.hvx2D(data.curVx,:),'Color','y','LineWidth',2);
end;

if (clicked3D),
    clickfig = data.fig3D;
else
    clickfig = data.fig2D;
end;

set(data.hNameEdit,'Visible','off');

%sort out double clicks or right clicks
c = get(clickfig,'CurrentPoint');
seltype = get(clickfig,'SelectionType');

switch seltype,
 case 'open',                           % double click
  set(data.hNameEdit,'Parent',clickfig,'Position',[c 100 20], ...
                    'String',data.vxr(ind).name, 'Visible','on');
  %disable the keypressfcn until they're done editing, because I often
  %open the edit, but forget to click in it and start typing immediately
  %and get weird effects
  %there doesn't seem to be any way just to set the focus directly to the
  %edit box, which would be the better solution here...
  set(clickfig,'KeyPressFcn','');

 case 'alt',                         % right click
  data = ucShowContextMenu(data,clickfig,c, ind);
end;
        
guidata(data.fig3D,data);

% -----------------------------------------------------------
function data = ucShowContextMenu(data,fig,pos, vxNum)

guidata(data.fig3D,data);
%set important parameters first
set(data.hVortexMenu,'Parent',fig, 'Position',pos, 'UserData',vxNum);
%then make it visible - this way we ensure that UserData really does
%equal vxNum before the menu comes up
set(data.hVortexMenu,'Visible','on');

%we'll step back in here after the menu callback, so update the data
%structure
data = guidata(data.fig3D);

% -----------------------------------------------------------
function ucMenuCallback(obj,eventdata,fcn)

h = guidata(obj);
if (ishandle(h)),
    data = guidata(h);
else
    data = h;
end;

vxNum = get(data.hVortexMenu,'UserData');

ind = find(data.vxr(vxNum).frames == data.curFrame);
if (isempty(ind)),                            % shouldn't ever happen
    return;
end;

vxrng = [];
switch fcn,
 case 'start',
  vxrng = [ind length(data.vxr(vxNum).frames)];

 case 'end',
  vxrng = [1 ind];

 case 'split',
  vxrng = [1 ind; ind+1 length(data.vxr(vxNum).frames)];

  %create a copy of the current vortex
  vx = data.vxr(vxNum);
  nvx = length(data.vxr) + 1;
  data.vxr(nvx) = vx;
  data.hvx3D(nvx) = copyobj(data.hvx3D(vxNum), ...
                            get(data.hvx3D(vxNum),'Parent'));
  data.hvx2D(nvx,:) = -1;
  data.vxrGood(nvx) = true;
  data.hvxlabel(nvx) = -1;

  %give it a new number and make sure it's not highlighted
  set(data.hvx3D(nvx),'UserData',nvx, 'EdgeColor','none');
  
  vxNum = [vxNum nvx];

 case 'delete',
  if (ishandle(data.hvx3D(data.curVx))),
      delete(data.hvx3D(data.curVx));
  end;
  if (ishandle(data.hvx2D(data.curVx))),
      delete(data.hvx2D(data.curVx,:));
  end;
  data.vxrGood(data.curVx) = false;
  data.curVx = [];

  data = ucUpdate(data, 'vortices');
end;

if (~isempty(vxrng)),
    for i = 1:size(vxrng,1),
        xd = get(data.hvx3D(vxNum(i)),'XData');
        yd = get(data.hvx3D(vxNum(i)),'YData');
        zd = get(data.hvx3D(vxNum(i)),'ZData');
        cd = get(data.hvx3D(vxNum(i)),'CData');

        vxrng1 = vxrng(i,1):vxrng(i,2);
        set(data.hvx3D(vxNum(i)),'XData',xd(:,vxrng1),...
                          'YData',yd(:,vxrng1),...
                          'ZData',zd(:,vxrng1),...
                          'CData',cd(:,vxrng1));

        %pull out the appropriate frames from all the fields in vxr
        len = length(data.vxr(vxNum(i)).frames);
        C = struct2cell(data.vxr(vxNum(i)));
        for j = 1:length(C),
            if (isnumeric(C{j}) & (size(C{j},2) == len)),
                C{j} = C{j}(vxrng1);
            end;
        end;
        S = cell2struct(C,fieldnames(data.vxr(vxNum(i))),1);

        data.vxr(vxNum(i)) = S;
    end;
end;

guidata(data.fig3D,data);

% -----------------------------------------------------------
function ucNameEditCallback(obj,eventdata)

h = guidata(obj);
if (ishandle(h)),
    data = guidata(h);
else
    data = h;
end;

data.vxr(data.curVx).name = get(obj,'String');
set(obj,'Visible','off');

if (ishandle(data.hvxlabel(data.curVx))),
    %change the existing label
    set(data.hvxlabel(data.curVx),...
        'String',data.vxr(data.curVx).name);
else
    %create a new one
    figure(data.fig2D);
    xd = get(data.hvx2D(data.curVx,1),'XData');
    yd = get(data.hvx2D(data.curVx,1),'YData');
    data.hvxlabel(data.curVx) = text(xd(end),yd(end), ...
                                     data.vxr(data.curVx).name,...
                                     'Color',data.vecCol,...
                                     'ButtonDownFcn',@ucSelectSurface);
end;
%change the transparency to indicate that the vortex is named
alpha(data.hvx3D(data.curVx),1);

%reset the keypressfcns
set(data.fig2D,'KeyPressFcn',@ucKeyPressFcn);
set(data.fig3D,'KeyPressFcn',@ucKeyPressFcn);

guidata(data.fig3D,data);

% -----------------------------------------------------------
function ucQuit(data)

% remove the various callbacks
set([data.fig3D data.fig2D],'KeyPressFcn','', ...
                  'CloseRequestFcn','closereq');
set(data.hvx2D(ishandle(data.hvx2D)), 'ButtonDownFcn','');
set(data.hvx3D(ishandle(data.hvx3D)), 'ButtonDownFcn','');

if (ishandle(data.hvx2D(data.curVx,1))),
    set(data.hvx2D(data.curVx,1),'Color','r','LineWidth',0.5);
end;
set(data.hvx3D(data.curVx),'EdgeColor','none');

delete(data.hSlice);
delete(data.hvx2D(ishandle(data.hvx2D(:,2)),2));

if (isempty(data.varnames{1})),
    data.varnames{1} = 'vxr';
end;
if (isempty(data.varnames{2})),
    data.varnames{2} = 'vxrGood';
end;
names = inputdlg({'Vortex structure:', 'Good vortex boolean vector:'},...
                 'Save variables?',1,data.varnames);
if (~isempty(names)),
    %sort the vortex structure, in case we split any vortices, because
    %they'll end up at the end of the vortex list, but should be in the
    %middle somewhere.
    %sort by start frame (ascending) and number of frames (descending)
    for i = 1:length(data.vxr),
        sortvals(i,1) = data.vxr(i).frames(1);
        sortvals(i,2) = -length(data.vxr(i).frames);
    end;
    [q,ord] = sortrows(sortvals);

    assignin('base',names{1},data.vxr(ord));
    assignin('base',names{2},data.vxrGood(ord));
end;


% -----------------------------------------------------------
function ucCloseFig(obj,eventdata)

h = guidata(obj);
if (ishandle(h)),
    data = guidata(h);
else
    data = h;
end;

ucQuit(data);

delete(obj);               

% -----------------------------------------------------------
function ind = findvxrinframe(vxr,fr)

vxind = 1;
ind = [];
while (vxind <= length(vxr)),
    if (any(vxr(vxind).frames == fr)),
        %mark vortices on this frame
        ind(end+1) = vxind;
    end;
    vxind = vxind+1;
end;



