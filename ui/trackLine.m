function [Xfr,Yfr] = trackLine(files,xfr,yfr, varargin)

opt.savecallback = [];
opt = parsevarargin(opt,varargin,4, 'typecheck',false);

if ~isempty(opt.savecallback) && ~isa(opt.savecallback, 'function_handle')
    error('savecallback must be a function handle');
end

if (nargin == 1),
  xfr = [];
  yfr = [];
end;

if (iscellstr(files)),
    handles.imType = 'cellstr';
    handles.imNames = files;
    handles.imNum = 1;
    handles.imTotal = length(handles.imNames);

    info = imfinfo(handles.imNames{1});
    handles.imWidth = info.Width;
    handles.imHeight = info.Height;

    handles.imFrameSkip = 1;
elseif (ischar(files)),
  info = VideoReader(files);
  handles.imType = 'avi';
  handles.imNames = files;
  handles.imNum = 1;
  handles.imTotal = info.Duration*info.FrameRate;
  handles.imWidth = info.Width;
  handles.imHeight = info.Height;
  handles.FrameRate = info.FrameRate;
  
  skip = inputdlg('Frame skip?','trackLine',1,{'1'});
  handles.imFrameSkip = str2num(skip{1});
else
    error('Unknown type of image sequence.');
end;

handles.savecallback = opt.savecallback;

handles.Figure = figure;
handles.Axes = axes('Parent',handles.Figure);

handles = createMenus(handles);

switch handles.imType,
 case 'cellstr',
  nm = handles.imNames{1};
 case 'avi',
  nm = handles.imNames;
end;
set(handles.Figure, 'Name', ...
                  sprintf('trackLine: %s (%d/%d)', nm,...
                          1,handles.imTotal),...
                  'KeyPressFcn', @KeyPress_Callback,...
                  'DoubleBuffer','on');

if (~isempty(xfr)),
    if (iscell(xfr) && (length(xfr) == handles.imTotal) && ...
        ~iscell(xfr{1})),
        for i = 1:length(xfr),
            xfr{i} = {xfr{i}};
            yfr{i} = {yfr{i}};
        end;
    elseif (isnumeric(xfr) && (size(xfr,2) == handles.imTotal)),
        xnum = xfr;
        ynum = yfr;
        xfr = cell(1,handles.imTotal);
        yfr = cell(1,handles.imTotal);
        for i = 1:size(xnum,2),
            xfr{i} = {xnum(:,i)};
            yfr{i} = {ynum(:,i)};
        end;
    end;

    handles.Xfr = xfr;
    handles.Yfr = yfr;

    handles.X = handles.Xfr{1}(1:end-1);
    handles.Y = handles.Yfr{1}(1:end-1);
    handles.Xcur = handles.Xfr{1}{end};
    handles.Ycur = handles.Yfr{1}{end};
else
  handles.Xcur = [];
  handles.Ycur = [];
  handles.X = {};
  handles.Y = {};
  handles.Xfr = cell(1,handles.imTotal);
  handles.Yfr = cell(1,handles.imTotal);
end;

handles.LineHandles = [];
handles.curLineHandle = -1;

handles.im = loadImage(handles);
handles = Update(handles,'image');

% Update handles structure
guidata(handles.Figure, handles);

% UIWAIT makes trackLine wait for user response (see UIRESUME)
uiwait(handles.Figure);

if (ishandle(handles.Figure)),
    handles = guidata(handles.Figure);

    Xfr = handles.Xfr;
    Yfr = handles.Yfr;
    Xfr{handles.imNum} = {handles.X{:} handles.Xcur};
    Yfr{handles.imNum} = {handles.Y{:} handles.Ycur};
    
    delete(handles.Figure);
else
    Xfr = [];
    Yfr = [];
end;

% --------------------------------------------------------------------
function handles = createMenus(handles)

imMenu = uicontextmenu;
uimenu(imMenu, 'Label', 'New line', 'Callback', @NewLine_Callback);

lnMenu = uicontextmenu;
uimenu(lnMenu, 'Label', 'Add point', 'Callback', @AddPoint_Callback);
uimenu(lnMenu, 'Label', 'Delete point', 'Callback', @DeletePoint_Callback);
uimenu(lnMenu, 'Label', 'New line', 'Callback', @NewLine_Callback, ...
       'Separator','on');
uimenu(lnMenu, 'Label', 'Delete line', 'Callback', @DeleteLine_Callback);

handles.ImageMenu = imMenu;
handles.LineMenu = lnMenu;

% --------------------------------------------------------------------
function handles = changeImage(handles, delta, fr)

if ~isempty(handles.Xcur)
    saveData(handles);
end

if (isfinite(delta)),
    i0 = handles.imNum;
    handles.imNum = handles.imNum + delta;
    i1 = handles.imNum;
else
    i0 = handles.imNum;
    handles.imNum = fr;
    i1 = fr;
end;

if ((handles.imNum < 1) || (handles.imNum > handles.imTotal)),
    beep;
    handles.imNum = i0;
else
    if (~isempty(handles.Xcur)),
        handles.Xfr{i0} = {handles.X{:} handles.Xcur};
        handles.Yfr{i0} = {handles.Y{:} handles.Ycur};
    end;

    if (~isempty(handles.Xfr{i1})),
        handles.X = handles.Xfr{i1}(1:end-1);
        handles.Y = handles.Yfr{i1}(1:end-1);
        handles.Xcur = handles.Xfr{i1}{end};
        handles.Ycur = handles.Yfr{i1}{end};

        if (ishandle(handles.curLineHandle)),
            delete(handles.curLineHandle);
        end;
        handles.curLineHandle = -1;
    end;

    handles.im = loadImage(handles);
    handles = Update(handles, 'image');
end;    


% --------------------------------------------------------------------
function saveData(handles)

if isa(handles.savecallback, 'function_handle')
    switch handles.imType
        case 'cellstr'
            imname = handles.imNames{handles.imNum};
        case 'avi'
            imname = handles.imNames;
    end
    feval(handles.savecallback, handles.Xcur, handles.Ycur, ...
        imname, handles.imNum);
end

% --------------------------------------------------------------------
function im = loadImage(imdata)

switch (imdata.imType),
 case 'cellstr',
  im = imread(imdata.imNames{imdata.imNum});
 case 'avi',
  vid = VideoReader(imdata.imNames);
  vid.CurrentTime = imdata.imNum/imdata.FrameRate;
  im = readFrame(vid);
end;

% --------------------------------------------------------------------
function handles = Update(handles, what)

if ischar(what)
    what = {what};
end

if ismember('image',what)
  if (isfield(handles,'imHandle') && ishandle(handles.imHandle)),
    set(handles.imHandle, 'CData', handles.im);
  else
    hold on;
    handles.imHandle = imshow(handles.im,'InitialMagnification','fit');
    set(handles.imHandle, 'ButtonDownFcn', @ImButtonDown_Callback,...
                      'UIContextMenu',handles.ImageMenu);
    hold off;
    axis equal ij on;
  end;

  switch handles.imType,
   case 'cellstr',
    nm = handles.imNames{1};
   case 'avi',
    nm = handles.imNames;
  end;
  set(handles.Figure, 'Name', ...
                    sprintf('%s %d/%d',nm,handles.imNum,handles.imTotal));

  % make sure the image is behind everything else
  cc = get(handles.Axes, 'Children');
  k = find(cc == handles.imHandle);
  cc = cc([1:k-1 k+1:end k]);
  set(handles.Axes, 'Children', cc);
  
  what = [what, 'curline','lines'];
end;

if ismember('curline',what)
    if (~isempty(handles.Xcur)),
        if (~isfield(handles,'curLineHandle') || ...
            ~ishandle(handles.curLineHandle)),
            handles.curLineHandle = addplot(handles.Xcur,handles.Ycur, ...
                                            'ro-');
            set(handles.curLineHandle, 'UIContextMenu', handles.LineMenu, ...
                'ButtonDownFcn',@LnButtonDown_Callback);
        else
            set(handles.curLineHandle, 'XData',handles.Xcur, ...
                              'YData',handles.Ycur);
        end;
    end;
    what = [what, 'lines'];
end;

if ismember('lines',what)
    if (~isempty(handles.X)),
        delete(handles.LineHandles);

        handles.LineHandles = [];
        for i = 1:length(handles.X),
            handles.LineHandles(i) = addplot(handles.X{i},handles.Y{i}, 'y.-');
        end;
        set(handles.LineHandles, 'UIContextMenu', handles.LineMenu, ...
             'ButtonDownFcn', @LnButtonDown_Callback);
    end;
end;

guidata(handles.Figure, handles);

% --------------------------------------------------------------------
function scale = getPtScale(ax)

u = get(ax, 'Units');
set(ax, 'Units', 'points');
pos = get(ax, 'Position');
set(ax, 'Units', u);

xl = get(ax, 'XLim');
yl = get(ax, 'YLim');

scx = pos(3)/diff(xl);
scy = pos(4)/diff(yl);

scale = max(scx,scy);


% --------------------------------------------------------------------
function [pt,nextpt,dist] = findSelPt(xs,ys, x,y, seldist)

dist = (xs - x).^2 + (ys - y).^2;
[dist,pts] = sort(dist);

if (nargin == 5),
    pt = pts(dist <= seldist^2);
    nextpt = [];
    dist = [];
else
    pt = pts(1);
    nextpt = pts(2);
    dist = sqrt(dist(1));
end;

% --------------------------------------------------------------------
function ImButtonDown_Callback(hObject, eventdata)

handles = guidata(hObject);

c = get(handles.Axes, 'CurrentPoint');
tp = get(handles.Figure, 'SelectionType');

switch (tp),
 case 'normal',
  handles.Xcur(end+1) = c(1,1);
  handles.Ycur(end+1) = c(1,2);
end;

guidata(hObject, handles);
Update(handles, 'curline');

% --------------------------------------------------------------------
function LnButtonDown_Callback(hObject, eventdata)

handles = guidata(hObject);

c = get(handles.Axes, 'CurrentPoint');
tp = get(handles.Figure, 'SelectionType');

switch (tp),
 case 'normal',
  % start dragging -- may not actually happen, if the button up comes 
  % quickly enough 
  ptScale = getPtScale(handles.Axes);
  sz = get(hObject, 'MarkerSize');
  pt = findSelPt(c(1,1),c(1,2), handles.Xcur,handles.Ycur, ...
                         1.1*sz/ptScale);
  if (~isempty(pt)),
      handles.dragType = 'point';
      handles.dragPoint = pt;
  else
      handles.dragType = 'line';
  end;

  % change line selection
  if (handles.curLineHandle ~= hObject),
      k = find(handles.LineHandles == hObject);
      if (isempty(k) || (length(k) > 1)),
          error('Weirdness!');
      end;

      set(handles.curLineHandle, 'Color','y', 'Marker','.');
      set(hObject, 'Color','r', 'Marker','o');

      handles.LineHandles = handles.LineHandles([1:k-1 k+1:end]);
      handles.LineHandles(end+1) = handles.curLineHandle;
      handles.curLineHandle = hObject;

      cx = handles.X{k};
      cy = handles.Y{k};
      handles.X = handles.X([1:k-1 k+1:end]);
      handles.Y = handles.Y([1:k-1 k+1:end]);
      handles.X{end+1} = handles.Xcur;
      handles.Y{end+1} = handles.Ycur;
      handles.Xcur = cx;
      handles.Ycur = cy;
  end;

  handles.dragStart = c(1,1:2);
  handles.Xcur0 = handles.Xcur;
  handles.Ycur0 = handles.Ycur;

  Update(handles,1);

  set(handles.Figure, 'WindowButtonMotionFcn',@MouseMove_Callback,...
         'WindowButtonUpFcn', @ButtonUp_Callback);
  guidata(handles.Figure,handles);
end;

% --------------------------------------------------------------------
function MouseMove_Callback(hObject, eventdata)

handles = guidata(hObject);

c = get(handles.Axes, 'CurrentPoint');

switch (handles.dragType),
 case 'point',
  handles.Xcur(handles.dragPoint) = c(1,1);
  handles.Ycur(handles.dragPoint) = c(1,2);
 case 'line',
  handles.Xcur = handles.Xcur0 + c(1,1) - handles.dragStart(1);
  handles.Ycur = handles.Ycur0 + c(1,2) - handles.dragStart(2);
end;

set(handles.curLineHandle, 'XData',handles.Xcur, 'YData',handles.Ycur);
guidata(hObject, handles);
drawnow;

% --------------------------------------------------------------------
function ButtonUp_Callback(hObject, eventdata)

handles = guidata(hObject);

c = get(handles.Axes, 'CurrentPoint');
c = c(1,1:2);

scale = getPtScale(handles.Axes);
if (sum((c - handles.dragStart).^2 > 1/scale)), % different point
    MouseMove_Callback(hObject, eventdata);
end;

set(handles.Figure, 'WindowButtonUpFcn','', 'WindowButtonMotionFcn','');

guidata(hObject, handles);

% --------------------------------------------------------------------
function KeyPress_Callback(hObject, eventdata)

handles = guidata(hObject);

c = get(handles.Figure, 'CurrentCharacter');

switch lower(c),
    case char(28),                         % left arrow
        handles = changeImage(handles, -handles.imFrameSkip);
    case {char(29),char(13)}               % right arrow, enter
        handles = changeImage(handles, +handles.imFrameSkip);
    case 'n',                              % new line
        NewLine_Callback(hObject, eventdata);
        handles = guidata(hObject);
    case 'g',                              % go to
        fr = inputdlg('Frame?','Go to frame');
        fr = str2num(fr{1});
        handles = changeImage(handles, NaN, fr);
    case 'c',                              % clear frame
        delete(handles.LineHandles);
        delete(handles.curLineHandle);
        handles.Xcur = [];
        handles.Ycur = [];
        handles.X = {};
        handles.Y = {};
        handles.Xfr{handles.imNum} = [];
        handles.Yfr{handles.imNum} = [];
    case 'z',                              % zoom
        zoom toggle;
    case 's',
        saveData(handles);
        
    case 'q',                              % quit
        saveData(handles);
        uiresume;
end;

guidata(hObject, handles);
  
% --------------------------------------------------------------------
function AddPoint_Callback(hObject, eventdata)

handles = guidata(hObject);

c = get(handles.Axes, 'CurrentPoint');

[pt1,pt2,dist] = findSelPt(c(1,1),c(1,2), handles.Xcur,handles.Ycur);
pt = [pt1 pt2];

% linearly interpolate a new point
totalDist = sqrt(diff(handles.Xcur(pt)).^2 + diff(handles.Ycur(pt)).^2);
ptx = handles.Xcur(pt(1)) + dist*diff(handles.Xcur(pt))/totalDist;
pty = handles.Ycur(pt(1)) + dist*diff(handles.Ycur(pt))/totalDist;

pt = sort(pt);
handles.Xcur = [handles.Xcur(1:pt(1)) ptx handles.Xcur(pt(2):end)];
handles.Ycur = [handles.Ycur(1:pt(1)) pty handles.Ycur(pt(2):end)];

guidata(hObject, handles);
Update(handles, 'curline');

% --------------------------------------------------------------------
function DeletePoint_Callback(hObject, eventdata)

handles = guidata(hObject);

c = get(handles.Axes, 'CurrentPoint');

ptScale = getPtScale(handles.Axes);
sz = get(handles.curLineHandle, 'MarkerSize');
pt = findSelPt(c(1,1),c(1,2), handles.Xcur,handles.Ycur, ...
               1.1*sz/ptScale);

if (~isempty(pt)),
    pt = pt(1);
    handles.Xcur = handles.Xcur([1:pt-1 pt+1:end]);
    handles.Ycur = handles.Ycur([1:pt-1 pt+1:end]);

    guidata(hObject, handles);
    Update(handles, 'curline');
else
    beep;
end;

% --------------------------------------------------------------------
function NewLine_Callback(hObject, eventdata)

handles = guidata(hObject);

handles.LineHandles(end+1) = handles.curLineHandle;
handles.curLineHandle = -1;

handles.X{end+1} = handles.Xcur;
handles.Y{end+1} = handles.Ycur;
handles.Xcur = [];
handles.Ycur = [];

guidata(hObject, handles);
Update(handles, 'lines');

% --------------------------------------------------------------------
function DeleteLine_Callback(hObject, eventdata)

handles = guidata(hObject);

delete(handles.curLineHandle);
handles.curLineHandle = -1;

handles.Xcur = handles.X{end};
handles.Ycur = handles.Y{end};
handles.X = handles.X(1:end-1);
handles.Y = handles.Y(1:end-1);

% we might have deleted the current object
guidata(handles.Figure, handles);
Update(handles, 'lines');





