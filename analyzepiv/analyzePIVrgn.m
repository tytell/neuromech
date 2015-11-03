function [data,rgnfcns] = analyzePIVrgn(data)

rgnfcns.apDrawRgn = @apDrawRgn;
rgnfcns.apDrawRgnButtDown = @apDrawRgnButtDown;
rgnfcns.apDrawRgnDrag = @apDrawRgnDrag;
rgnfcns.apDrawRgnButtUp = @apDrawRgnButtUp;
rgnfcns.apDrawRgnLine = @apDrawRgnLine;
rgnfcns.apDrawRgnSquare = @apDrawRgnSquare;
rgnfcns.apDrawRgnCircle = @apDrawRgnCircle;
rgnfcns.apDrawPolyRgn = @apDrawPolyRgn;
rgnfcns.apDrawPolyRgnDrag = @apDrawPolyRgnDrag;
rgnfcns.apDrawPolyRgnClick = @apDrawPolyRgnClick;
rgnfcns.apClickRgn = @apClickRgn;
rgnfcns.apDragRgnMove = @apDragRgnMove;
rgnfcns.apClickRgnButtUp = @apClickRgnButtUp;
rgnfcns.apSetCurRgn = @apSetCurRgn;
rgnfcns.apAddRgn = @apAddRgn;
rgnfcns.apStopRgn = @apStopRgn;
rgnfcns.apDeleteRgn = @apDeleteRgn;
rgnfcns.apSetRgnFrame = @apSetRgnFrame;
rgnfcns.apInterpRgn = @apInterpRgn;
rgnfcns.apIDVortices = @apIDVorticesCallback;
rgnfcns.apLoadRgns = @apLoadRgns;

data.rgnNameNum = 0;

set(data.RgnLineButton,'Callback',{@apDrawRgn,'line'});
set(data.RgnSquareButton,'Callback',{@apDrawRgn,'square'});
set(data.RgnCircleButton,'Callback',{@apDrawRgn,'circle'});
set(data.RgnPolygonButton,'Callback',@apDrawPolyRgn);
set(data.IDVorticesButton,'Callback',@apIDVorticesCallback);

set(data.RgnMenu,'Callback',{@apSetupRgnMenu,data.Panel}, ...
                 'Parent',data.Figure);
set(data.KeyframeRgnMenu,'Callback',{@apKeyframeRgnMenu,data.Panel});
set(data.SetRgnStartMenu,'Callback',{@apSetRgnStartMenu,data.Panel});
set(data.SetRgnStopMenu,'Callback',{@apSetRgnStopMenu,data.Panel});
set(data.DeleteRgnMenu,'Callback',{@apDeleteRgnMenu,data.Panel});

set(data.LoadRgnButton,'Callback',@apLoadRgns);

% -------------------------------------------------
function apLoadRgns(obj, eventdata)

data = guidata(obj);

if (size(data.Regions,1) > 0),
    switch questdlg('Overwrite current regions or append?',...
                    'Current regions','Yes','Append','Cancel','Append'),
     case 'Yes',
      data.Regions = data.Regions(1:0,:);
     case 'Append',
      % do nothing for the moment
     case 'Cancel',
      return;
    end;
end;

[file,path] = uigetfile({'*.mat','Matlab data file (*.mat)'},...
                                'Load regions from file...');

file = fullfile(path,file);
load(file);

if (exist('Regions')),
    fprintf('Found %d regions and %d frames.\n',size(Regions));
    if (size(Regions,2) ~= data.nFrames),
        warndlg('Cannot import regions with different numbers of frames.');
        return;
    end;

    axes(data.Axes);
    newrgn1 = size(data.Regions,1);
    for i = 1:size(Regions,1),
        switch Regions(i,1).type,
         case 'line',
          tp = 'Ln';
         case 'circle',
          tp = 'Cl';
         case 'square',
          tp = 'Sq';
        end;
        Regions(i,data.curFrame).handle = ...
            line('XData',Regions(i,data.curFrame).x,...
                 'YData',Regions(i,data.curFrame).y,...
                 'Color','k',...
                 'ButtonDownFcn', {@apClickRgn, data.Panel},...
                 'Tag',tp); 

        data.rgnNames{i+newrgn1} = sprintf('%s %d',tp,i+newrgn1);
        Regions(i,data.curFrame).label = ...
            text(Regions(i,data.curFrame).x(end),...
                 Regions(i,data.curFrame).y(end),...
                 data.rgnNames{i+newrgn1},...
                 'HorizontalAlignment','left',...
                 'Color','r','HitTest','off');
    end;
    data.Regions(end+(1:size(Regions,1)),:) = Regions;
    data.rgnNameNum = size(data.Regions,1);

    if (~isfield(data,'curRgn') | isnan(data.curRgn)),
        data = apSetCurRgn(data,newrgn1+1);
    end;
    data = apSetRgnFrame(data, data.curFrame,data.curFrame);
    data = analyzePIVupdate(data, 'regions',{'recalc','rgn',data.curRgn});
end;

guidata(obj,data);

% -------------------------------------------------
function apIDVorticesCallback(obj, eventdata)

data = guidata(obj);
data = analyzePIVvortices(data);
guidata(obj,data);


% -------------------------------------------------
function apDrawRgn(obj, eventdata, rgntype)

data = guidata(obj);

if (get(obj,'Value')),
    switch rgntype,
     case 'line',
      drawFcn = @apDrawRgnLine;
      tp = 'Ln';
      set([data.RgnSquareButton data.RgnCircleButton ...
           data.RgnPolygonButton],'Value',0);
     case 'circle',
      drawFcn = @apDrawRgnCircle;
      tp = 'Cl';
      set([data.RgnLineButton data.RgnSquareButton ...
           data.RgnPolygonButton],'Value',0);
     case 'square',
      drawFcn = @apDrawRgnSquare;
      tp = 'Sq';
      set([data.RgnLineButton data.RgnCircleButton ...
           data.RgnPolygonButton],'Value',0);
    end;

    axes(data.Axes);
    data.dragHandle = line('XData',[],'YData',[],...
                           'Color','k',...
                           'Tag',tp); 
% ,'UIContextMenu',data.RgnMenu);
    data.rgnTool = obj;
    data.rgnType = rgntype;

    set(data.Figure,'WindowButtonDownFcn',...
                    {@apDrawRgnButtDown,drawFcn,data.Panel});
    if (~isempty(data.Regions)),
        set(cat(1,data.Regions(:,data.curFrame).handle),'HitTest','off');
        if (ishandle(data.rgnSzHandle)),
            set([data.rgnSzHandle data.rgnRotHandle],'HitTest','off');
        end;
    end;
else
    if (ishandle(data.dragHandle)),
        delete(data.dragHandle);
    end;
    set(data.Figure,'WindowButtonDownFcn','');
    if (~isempty(data.Regions)),
        set(cat(1,data.Regions(:,data.curFrame).handle),'HitTest','on');
        if (ishandle(data.rgnSzHandle)),
            set([data.rgnSzHandle data.rgnRotHandle],'HitTest','on');
        end;
    end;
end;
guidata(obj,data);


% -------------------------------------------------
function apDrawRgnButtDown(obj, eventdata, drawFcn, panel)

data = guidata(panel);

pt = get(data.Axes,'CurrentPoint');
data.dragStart = pt(1,1:2);

data.dragHandle = feval(drawFcn, data.dragHandle,data.dragStart,...
                        pt(1,1:2),'b',0);
data.OldMotionFcn = get(data.Figure,'WindowButtonMotionFcn');
set(data.Figure,'WindowButtonMotionFcn',...
                {@apDrawRgnDrag,drawFcn,data.Panel}, ...
                'WindowButtonUpFcn', ...
                {@apDrawRgnButtUp,drawFcn,data.Panel});

guidata(panel,data);

% -------------------------------------------------
function apDrawRgnDrag(obj, eventdata, drawFcn, panel)

data = guidata(panel);

pt = get(data.Axes,'CurrentPoint');
c = get(obj,'SelectionType');
data.dragHandle = feval(drawFcn, data.dragHandle,data.dragStart,...
                        pt(1,1:2),'b',strcmp(c,'alt'));

guidata(panel,data);

% -------------------------------------------------
function apDrawRgnButtUp(obj, eventdata, drawFcn, panel)

data = guidata(panel);

pt = get(data.Axes,'CurrentPoint');
c = get(obj,'SelectionType');
h = feval(drawFcn, data.dragHandle,data.dragStart,pt(1,1:2),...
      'r',strcmp(c,'alt'));
set(h, 'ButtonDownFcn', {@apClickRgn, panel});

data.dragHandle = -1;

set(data.Figure,'WindowButtonDownFcn','',...
                'WindowButtonMotionFcn',data.OldMotionFcn,...
                'WindowButtonUpFcn', '');
set(data.rgnTool, 'Value', 0);

data = apAddRgn(data, h);

set(cat(1,data.Regions(:,data.curFrame).handle),'HitTest','on');
if (ishandle(data.rgnSzHandle)),
    set([data.rgnSzHandle data.rgnRotHandle],'HitTest','on');
end;

guidata(panel,data);

% -------------------------------------------------
function h = apDrawRgnLine(h, pos1,pos2, col, constrain)

set(h,'XData',[pos1(1) pos2(1)], 'YData',[pos1(2) pos2(2)], ...
      'Color',col);

% -------------------------------------------------
function h = apDrawRgnSquare(h, pos1,pos2, col, constrain)

x = [pos1(1) pos2(1)];
y = [pos1(2) pos2(2)];

if (constrain),
    ax = min(abs(diff(x)),abs(diff(y)));
    x(2) = x(1) + ax*sign(diff(x));
    y(2) = y(1) + ax*sign(diff(y));
end;

x = x([1 2 2 1 1]);
y = y([1 1 2 2 1]);

set(h,'XData',x, 'YData',y, 'Color',col);


% -------------------------------------------------
function h = apDrawRgnCircle(h, pos1,pos2, col, constrain)

x0 =  [1.0000    0.9766    0.9076    0.7961    0.6474  ...
       0.4684    0.2675    0.0541   -0.1618   -0.3701  ...
       -0.5612   -0.7260   -0.8569   -0.9477   -0.9941 ...
       -0.9941   -0.9477   -0.8569   -0.7260   -0.5612 ...
       -0.3701   -0.1618    0.0541    0.2675    0.4684 ...
       0.6474    0.7961    0.9076    0.9766    1.0000];
y0 =  [0    0.2150    0.4199    0.6052    0.7622    ...
       0.8835    0.9635    0.9985    0.9868    0.9290  ...
       0.8277    0.6877    0.5156    0.3193    0.1081  ...
       -0.1081   -0.3193   -0.5156   -0.6877   -0.8277 ...
       -0.9290   -0.9868   -0.9985   -0.9635   -0.8835 ...
       -0.7622   -0.6052   -0.4199   -0.2150   0.0000];

ax1 = abs(pos2(1) - pos1(1));
ax2 = abs(pos2(2) - pos1(2));

if (constrain),
    ax1 = sqrt((pos2(1)-pos1(1)).^2 + (pos2(2)-pos1(2)).^2);
    ax2 = ax1;
end;

x = x0*ax1 + pos1(1);
y = y0*ax2 + pos1(2);
x(end) = x(1);
y(end) = y(1);

set(h,'XData',x, 'YData',y, 'Color',col);

% -------------------------------------------------
function apDrawPolyRgn(obj, eventdata)

if (get(obj,'Value')),
    data = guidata(obj);

    set([data.RgnLineButton data.RgnCircleButton ...
         data.RgnSquareButton],'Value',0);

    axes(data.Axes);
    data.dragHandle = line('XData',[],'YData',[],...
                           'Color','k',...
                           'Tag','Pg');
    data.rgnType = 'poly';

    set(data.Figure,'WindowButtonUpFcn', {@apDrawPolyRgnClick,data.Panel});
    if (~isempty(data.Regions)),
        set(cat(1,data.Regions(:,data.curFrame).handle),'HitTest','off');
        if (ishandle(data.rgnSzHandle)),
            set([data.rgnSzHandle data.rgnRotHandle],'HitTest','off');
        end;
    end;
else
    if (ishandle(data.dragHandle)),
        delete(data.dragHandle);
    end;
    set(data.Figure,'WindowButtonDownFcn','');
    if (~isempty(data.Regions)),
        set(cat(1,data.Regions(:,data.curFrame).handle),'HitTest','on');

        if (ishandle(data.rgnSzHandle)),
            set([data.rgnSzHandle data.rgnRotHandle],'HitTest','on');
        end;
    end;
end;

guidata(obj, data);

% -------------------------------------------------
function apDrawPolyRgnDrag(obj, eventdata, panel)

data = guidata(panel);

pt = get(data.Axes,'CurrentPoint');

xd = get(data.dragHandle,'XData');
yd = get(data.dragHandle,'YData');

if (length(xd) <= 2),
    xd(end) = pt(1,1);
    yd(end) = pt(1,2);
else
    xd(end-1) = pt(1,1);
    yd(end-1) = pt(1,2);
end;

set(data.dragHandle,'XData',xd,'YData',yd);

guidata(panel, data);

% -------------------------------------------------
function apDrawPolyRgnClick(obj, eventdata, panel)

data = guidata(panel);

pt = get(data.Axes,'CurrentPoint');
c = get(data.Figure,'SelectionType');

xd = get(data.dragHandle,'XData');
yd = get(data.dragHandle,'YData');

if (length(xd) == 0),                   % first click
    xd = pt(1,[1 1]);
    yd = pt(1,[2 2]);
    data.OldMotionFcn = get(data.Figure,'WindowButtonMotionFcn');
    set(data.Figure, 'WindowButtonMotionFcn', ...
                     {@apDrawPolyRgnDrag,data.Panel});
elseif (length(xd) == 2),
    xd = [xd(1:end) pt(1,1) xd(1)];
    yd = [yd(1:end) pt(1,2) yd(1)];
else
    xd = [xd(1:end-1) pt(1,1) xd(1)];
    yd = [yd(1:end-1) pt(1,2) yd(1)];
end;

if (~strcmp(c,'normal')),
    h = data.dragHandle;
    data.dragHandle = -1;

    % double click
    if ((xd(end-1) == xd(end-2)) & (yd(end-1) == yd(end-2))),
        xd = [xd(1:end-3) xd(1)];
        yd = [yd(1:end-3) yd(1)];
    end;
    set(h,'XData',xd,'YData',yd,'Color','r',...
          'ButtonDownFcn', {@apClickRgn, panel});

    data = apAddRgn(data, h);

    set(data.Figure,'WindowButtonUpFcn', '', ...
                    'WindowButtonMotionFcn', data.OldMotionFcn);
    set(data.RgnPolygonButton,'Value',0);

    set(cat(1,data.Regions(:,data.curFrame).handle),'HitTest','on');
    if (ishandle(data.rgnSzHandle)),
        set([data.rgnSzHandle data.rgnRotHandle],'HitTest','on');
    end;
else
    set(data.dragHandle,'XData',xd,'YData',yd);
end;
    
guidata(panel, data);



% -------------------------------------------------
function apClickRgn(obj, eventdata, panel)

data = guidata(panel);

switch get(data.Figure,'SelectionType'),
 case 'normal',
  c = get(data.Axes, 'CurrentPoint');
  data.dragStart = c(1,1:2);
  data.dragHandle = obj;
  data.dragX0 = get(obj,'XData');
  data.dragY0 = get(obj,'yData');

  data.OldMotionFcn = get(data.Figure,'WindowButtonMotionFcn');
  set(data.Figure,'WindowButtonMotionFcn',...
                  {@apDragRgnMove,panel}, ...
                  'WindowButtonUpFcn', ...
                  {@apClickRgnButtUp,panel});
 case 'alt',
  c = get(data.Figure,'CurrentPoint');
  data = apSetupRgnMenu(data,obj);
  
  % seems that we end up stepping back into here after we click
  % something in the popup menu
  % since we've almost always modified the data structure, we
  % need to be careful to update it

  guidata(panel,data);
  set(data.RgnMenu,'Position',c, 'Visible','on');
  data = guidata(panel);
end;

guidata(panel,data);

% -------------------------------------------------
function apDragRgnMove(obj, eventdata, panel)

data = guidata(panel);

c = get(data.Axes, 'CurrentPoint');

dx = c(1,1) - data.dragStart(1);
dy = c(1,2) - data.dragStart(2);

try
    set(data.dragHandle,'XData',data.dragX0+dx, 'YData',data.dragY0+dy);
catch err
    set(data.Figure,'WindowButtonMotionFcn','',...
                    'WindowButtonUpFcn', '');
    
    fprintf('Weird problem:\n');
    disp(err);
end 
    
guidata(panel,data);

% -------------------------------------------------
function apClickRgnButtUp(obj, eventdata, panel)

data = guidata(panel);

c = get(data.Axes, 'CurrentPoint');

dx = c(1,1) - data.dragStart(1);
dy = c(1,2) - data.dragStart(2);

set(data.dragHandle,'XData',data.dragX0+dx, 'YData',data.dragY0+dy);
set(data.Figure,'WindowButtonMotionFcn',data.OldMotionFcn,...
                'WindowButtonUpFcn', '');

data = apSetCurRgn(data, data.dragHandle);
data = analyzePIVupdate(data,'regions',{'recalc','rgn',data.curRgn});

guidata(panel,data);

% -------------------------------------------------
function apResizeRgn(obj, eventdata, panel)

data = guidata(panel);

pt = get(data.Axes,'CurrentPoint');
szptx = get(obj,'XData');
szpty = get(obj,'YData');

[q,ptnum0] = min((pt(1,1)-szptx).^2 + (pt(1,2)-szpty).^2);

ptnum = mod((ptnum0-1)+4,8)+1;           % get the opposite point
data.dragStart = [szptx(ptnum) szpty(ptnum)];
data.dragDir = sign([szptx(ptnum)-szptx(ptnum0) ...
                    szpty(ptnum)-szpty(ptnum0)]);

data.dragX0 = get(data.Regions(data.curRgn,data.curFrame).handle, ...
                  'XData');
data.dragY0 = get(data.Regions(data.curRgn,data.curFrame).handle, ...
                  'YData');

data.OldMotionFcn = get(data.Figure,'WindowButtonMotionFcn');
set(data.Figure,'WindowButtonMotionFcn',...
                {@apResizeRgnMove,panel}, ...
                'WindowButtonUpFcn', ...
                {@apResizeRgnButtUp,panel});

guidata(panel,data);

% -------------------------------------------------
function apResizeRgnMove(obj, eventdata, panel)

data = guidata(panel);

pt = get(data.Axes,'CurrentPoint');

sz = data.dragStart - pt(1,1:2);

szptx = get(data.rgnSzHandle,'XData');
szpty = get(data.rgnSzHandle,'YData');

sz = sz./[range(data.dragX0) range(data.dragY0)].*data.dragDir;
sz(sz == 0) = 1;

xd = (data.dragX0 - data.dragStart(1))*sz(1) + data.dragStart(1);
yd = (data.dragY0 - data.dragStart(2))*sz(2) + data.dragStart(2);

set(data.Regions(data.curRgn,data.curFrame).handle,...
    'XData',xd,'YData',yd);

xd = [min(xd) (min(xd)+max(xd))/2 max(xd)];
yd = [min(yd) (min(yd)+max(yd))/2 max(yd)];
szptx = xd([1 1 1 2 3 3 3 2]);
szpty = yd([1 2 3 3 3 2 1 1]);
set(data.rgnSzHandle,'XData',szptx,'YData',szpty);

guidata(panel,data);

% -------------------------------------------------
function apResizeRgnButtUp(obj, eventdata, panel)

data = guidata(panel);

set(data.Figure,'WindowButtonMotionFcn', data.OldMotionFcn, ...
                'WindowButtonUpFcn', '');

hand = data.Regions(data.curRgn,data.curFrame).handle;
data = apSetCurRgn(data, hand);
data = analyzePIVupdate(data, 'regions',{'recalc','rgn',data.curRgn});

guidata(panel,data);


% -------------------------------------------------
function apRotateRgn(obj, eventdata, panel)

data = guidata(panel);

pt = get(data.Axes,'CurrentPoint');

data.dragStart = pt(1,1:2);

data.dragX0 = get(data.Regions(data.curRgn,data.curFrame).handle, ...
                  'XData');
data.dragY0 = get(data.Regions(data.curRgn,data.curFrame).handle, ...
                  'YData');

data.OldMotionFcn = get(data.Figure,'WindowButtonMotionFcn');
set(data.Figure,'WindowButtonMotionFcn',...
                {@apRotateRgnMove,panel}, ...
                'WindowButtonUpFcn', ...
                {@apRotateRgnButtUp,panel});

guidata(panel,data);

% -------------------------------------------------
function apRotateRgnMove(obj, eventdata, panel)

data = guidata(panel);

pt = get(data.Axes,'CurrentPoint');

rgnctr = [(min(data.dragX0)+max(data.dragX0))/2 ...
          (min(data.dragY0)+max(data.dragY0))/2];

ang = atan2(pt(1,1)-rgnctr(1),rgnctr(2)-pt(1,2));

xd = cos(ang)*(data.dragX0-rgnctr(1)) - sin(ang)*(data.dragY0-rgnctr(2)) ...
     + rgnctr(1);
yd = sin(ang)*(data.dragX0-rgnctr(1)) + cos(ang)*(data.dragY0-rgnctr(2)) ...
     + rgnctr(2);

set(data.Regions(data.curRgn,data.curFrame).handle,...
    'XData',xd,'YData',yd);

guidata(panel,data);

% -------------------------------------------------
function apRotateRgnButtUp(obj, eventdata, panel)

data = guidata(panel);

set(data.Figure,'WindowButtonMotionFcn', data.OldMotionFcn, ...
                'WindowButtonUpFcn', '');

% TODO: small changes can cause mirroring, somehow, sometimes.  Fix
% this!
hand = data.Regions(data.curRgn,data.curFrame).handle;
data = apSetCurRgn(data, hand);
data = analyzePIVupdate(data, 'regions',{'recalc','rgn',data.curRgn});

guidata(panel,data);

% -------------------------------------------------
function data = apSetCurRgn(data, rgn)

rgns = data.Regions(:,data.curFrame);

if (~ishandle(rgn)),
    data.curRgn = NaN;
    if ishandle(data.rgnSzHandle)
        delete(data.rgnSzHandle);
    end
    if ishandle(data.rgnRotHandle)
        delete(data.rgnRotHandle);
    end
    return;
end;

if ((rgn >= 1) && (rgn <= length(rgns)) && ...
    (mod(rgn,1) == 0)),
    ind = rgn;
    updatefromhandle = 0;
elseif (ishandle(rgn)),
    ind = find(cat(1,rgns.handle) == rgn);
    updatefromhandle = 1;
end;

data.curRgn = ind;

rest = [1:ind-1 ind+1:length(rgns)];

rest = rest(find(cat(1,rgns(rest).status) == 2));
set(cat(1,rgns(rest).handle), 'Color', 'k');
set(cat(1,rgns(rest).label), 'Color', 'k');

if (updatefromhandle),
    xd = get(rgns(ind).handle,'XData');
    yd = get(rgns(ind).handle,'YData');
    rgns(ind).x = xd;
    rgns(ind).y = yd;
else
    xd = rgns(ind).x;
    yd = rgns(ind).y;
end;
rgns(ind).status = 2;

set(rgns(ind).handle, 'Color', 'r');
if (ishandle(rgns(ind).label)),
    x = rgns(ind).x(end);
    y = rgns(ind).y(end);
    set(rgns(ind).label, 'Color','r', ...
                      'Position',[x y 0]);
end;

data.Regions(:,data.curFrame) = rgns;

xd = [min(xd) (min(xd)+max(xd))/2 max(xd)];
yd = [min(yd) (min(yd)+max(yd))/2 max(yd)];
szptx = xd([1 1 1 2 3 3 3 2]);
szpty = yd([1 2 3 3 3 2 1 1]);

axes(data.Axes);
if (ishandle(data.rgnSzHandle)),
    set(data.rgnSzHandle,'XData',szptx,'YData',szpty);
else
    data.rgnSzHandle = line('XData',szptx,'YData',szpty,'Marker','s',...
                            'MarkerFaceColor','k',...
                            'MarkerEdgeColor','none',...
                            'LineStyle','none',...
                            'ButtonDownFcn',{@apResizeRgn,data.Panel});
end;

% figure out an absolute point displacement for the rotate point
set(data.Axes,'Units','points');
pos = get(data.Axes,'Position');
set(data.Axes,'Units','normalized');
ax = axis;
pointscale = (ax(2)-ax(1))/pos(3);      % axis units per point (in x)

rotptx = xd(2);
switch data.axisOrient,
 case 'ij',
  rotpty = yd(1)-10*pointscale;
 case 'xy',
  rotpty = yd(3)+10*pointscale;
end;

if (ishandle(data.rgnRotHandle)),
    set(data.rgnRotHandle,'XData',rotptx,'YData',rotpty);
else
    data.rgnRotHandle = line('XData',rotptx,'YData',rotpty,'Marker','s',...
                             'MarkerFaceColor','g',...
                             'MarkerEdgeColor','none',...
                             'LineStyle','none',...
                             'ButtonDownFcn',{@apRotateRgn,data.Panel});
end;

% -------------------------------------------------
function data = apAddRgn(data, h)

x1 = get(h,'XData');
y1 = get(h,'YData');
% check for size 0 regions (accidental quick button click)
if ((range(x1) == 0) | (range(y1) == 0)),
    return;
end;

cur = size(data.Regions,1)+1;
data.Regions(cur,1) = struct('handle',-1, 'type',data.rgnType, ...
                             'label',-1, 'x',[],'y',[],'status',0);
data.Regions(cur,2:data.nFrames) = ...
    repmat(data.Regions(cur,1),[1 data.nFrames-1]);

rgn = data.Regions(cur,data.curFrame);

rgn.handle = h;
rgn.x = x1;
rgn.y = y1;
rgn.type = data.rgnType;

rgn.status = 2;

% we keep a separate index for the number attached to a region name
% and we never decrease it, even when we delete a region, to make sure
% we don't accidentally create multiple regions with the same name
data.rgnNameNum = data.rgnNameNum + 1;

tp = get(h, 'Tag');
data.rgnNames{cur} = sprintf('%s %d',tp,data.rgnNameNum);
rgn.label = text(rgn.x(end),rgn.y(end),...
                 data.rgnNames{cur},...
                 'HorizontalAlignment','left',...
                 'Color','r','HitTest','off');
data.Regions(cur,data.curFrame) = rgn;

data.Regions(cur,data.curFrame+1).status = 1;
data.Regions(cur,data.curFrame+2:data.nFrames) = ...
    repmat(data.Regions(cur,data.curFrame+1),...
           [1 data.nFrames-data.curFrame-1]);

% add a row for this region in the calculation window
data = feval(data.calcGuiFcns.apAddCalc,data);

data = apSetCurRgn(data,cur);
data = analyzePIVupdate(data,'regions',{'recalc','rgn',cur});

% -------------------------------------------------
function data = apStopRgn(data,rgn)

for i = data.curFrame+1:data.nFrames,
    data.Regions(rgn,i).status = 0;
end;
data.Regions(rgn,data.curFrame).status = 2;

data = analyzePIVupdate(data,'regions',{'recalc','rgn',rgn});

% -------------------------------------------------
function data = apStartRgn(data,rgn)

for i = 1:data.curFrame-1,
    data.Regions(rgn,i).status = 0;
end;
data.Regions(rgn,data.curFrame).status = 2;

data = analyzePIVupdate(data,'regions',{'recalc','rgn',rgn});

% -------------------------------------------------
function data = apDeleteRgn(data,rgnnum)

% remove the region entirely
rgn = data.Regions(rgnnum,data.curFrame);
data.Regions = data.Regions([1:rgnnum-1 rgnnum+1:end],:);
data.rgnNames = data.rgnNames([1:rgnnum-1 rgnnum+1:end]);

% delete the graphics items
delete(rgn.handle);
delete(rgn.label);

% and reset the current region if this one was the current region
if (rgnnum == data.curRgn),
    data = apSetCurRgn(data,NaN);
end;

data = feval(data.calcGuiFcns.apDeleteCalc,data,rgnnum);
data = analyzePIVupdate(data,'regions',{'recalc'});

% -------------------------------------------------
function data = apSetupRgnMenu(data, rgnHandle)

% this function should only be called when the menu is popping up

rgn = find(cat(1,data.Regions(:,data.curFrame).handle) == rgnHandle);
switch (data.Regions(rgn,data.curFrame).status),
 case 2,
    set(data.KeyframeRgnMenu,'Checked','on');
 case 1,
    set(data.KeyframeRgnMenu,'Checked','off');
end;

% we need to keep track of which region we're operating on (it's not
% necessarily the current region), so we save it in UserData
set(data.RgnMenu,'UserData',rgn);


% -------------------------------------------------
function apKeyframeRgnMenu(obj, eventdata, panel)

data = guidata(panel);

% get the region we popped up for
rgn = get(data.RgnMenu, 'UserData');

chk = get(obj,'Checked');
if (chk),
    data.Regions(rgn,data.curFrame).status = 2;
else
    data.Regions(rgn,data.curFrame).status = 1;
end;
% use apSetRgnFrame to update the regions, not to change frames
apSetRgnFrame(data, data.curFrame,data.curFrame);

guidata(panel,data);

% -------------------------------------------------
function apSetRgnStartMenu(obj, eventdata, panel)

data = guidata(panel);

% get the region we popped up for
rgn = get(data.RgnMenu, 'UserData');
data = apStartRgn(data,rgn);

guidata(panel,data);

% -------------------------------------------------
function apSetRgnStopMenu(obj, eventdata, panel)

data = guidata(panel);

% get the region we popped up for
rgn = get(data.RgnMenu, 'UserData');
data = apStopRgn(data,rgn);

guidata(panel,data);


% -------------------------------------------------
function apDeleteRgnMenu(obj, eventdata, panel)

data = guidata(panel);

% get the region we popped up for
rgn = get(data.RgnMenu, 'UserData');
data = apDeleteRgn(data,rgn);

guidata(panel,data);

% -------------------------------------------------
function data = apSetRgnFrame(data, to,from)

data = apInterpRgn(data);

if (~isempty(data.Regions)),
    hand = cat(1,data.Regions(:,from).handle);
    label = cat(1,data.Regions(:,from).label);

    for i = 1:size(data.Regions,1),
        data.Regions(i,from).handle = -1;
        data.Regions(i,from).label = -1;

        data.Regions(i,to).handle = hand(i);
        data.Regions(i,to).label = label(i);
        if (data.Regions(i,to).status == 0),
            set(hand(i),'Visible','off');
            set(label(i),'Visible','off');
        else
            if (data.Regions(i,to).status == 1),
                col = [0.5 0.5 0.5];
                ls = '--';
                width = 0.5;
            else
                col = [0 0 0];
                ls = '-';
                width = 1;
            end;

            set(hand(i),'XData',data.Regions(i,to).x, ...
                'YData',data.Regions(i,to).y, 'Color',col,...
                'LineStyle',ls,'LineWidth',width,'Visible','on');
            set(label(i),'Position',...
                [data.Regions(i,to).x(end) data.Regions(i,to).y(end) ...
                0], 'Color',col,...
                'LineStyle',ls,'LineWidth',width,'Visible','on');
        end;
    end;

    if (~isnan(data.curRgn) & data.Regions(data.curRgn,to).status == 2),
        data = apSetCurRgn(data,data.curRgn);
    else
        k = find(cat(1,data.Regions(:,to).status) == 2);
        if (~isempty(k)),
            data = apSetCurRgn(data,k(end));
        else
            if (ishandle(data.rgnSzHandle)),
                delete(data.rgnSzHandle);
            end;
            if (ishandle(data.rgnRotHandle)),
                delete(data.rgnRotHandle);
            end;
        end;
    end;
end;    

% -------------------------------------------------
function data = apInterpRgn(data)

for i = 1:size(data.Regions,1),
    k = find(cat(2,data.Regions(i,:).status) == 2);
    kdel = find(cat(2,data.Regions(i,:).status) == 0);
    
    if (k(end) < data.nFrames),
        m = min(kdel(kdel > k(end)));
        if (isempty(m)),
            m = data.nFrames+1;
        end;
        m = m-1;

        if (m > k(end)),
            data.Regions(i,m).x = data.Regions(i,k(end)).x;
            data.Regions(i,m).y = data.Regions(i,k(end)).y;
            data.Regions(i,m).status = 1;
            k(end+1) = m;
        end;
    end;

    for j = 1:length(k)-1,
        xd1 = data.Regions(i,k(j)).x;
        yd1 = data.Regions(i,k(j)).y;
        xd2 = data.Regions(i,k(j+1)).x;
        yd2 = data.Regions(i,k(j+1)).y;

        n = k(j+1)-k(j);
        dx = (xd2-xd1)/n;
        dy = (yd2-yd1)/n;
        
        for a = 1:n-1,
            data.Regions(i,k(j)+a).x = xd1 + a*dx;
            data.Regions(i,k(j)+a).y = yd1 + a*dy;
        end;
    end;
end;

