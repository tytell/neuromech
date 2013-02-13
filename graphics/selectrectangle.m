function finalrect = selectrectangle(varargin)

opt.color = 'r';

if ((nargin >= 1) && (numel(varargin{1}) == 1) && ishandle(varargin{1}) && ...
        strcmp(get(varargin{1},'Type'), 'axes'))
    hax = varargin{1};
    p = 2;
else
    hax = gca;
    p = 1;
end;

initialrect = [];
if ((nargin >= p) && isnumeric(varargin{p}) && (numel(varargin{p}) == 4))
    initialrect = varargin{p};
    p = p+1;
end;

opt = parsevarargin(opt,varargin(p:end), p);

if (isempty(initialrect))
    axlim = axis(hax);
    
    w = axlim(2) - axlim(1);
    h = axlim(4) - axlim(3);
    
    initialrect = axlim + [w/4 -w/4 h/4 -h/4];
end;

w = initialrect(2) - initialrect(1);
h = initialrect(4) - initialrect(3);

xr = initialrect([1 2 2 1 1]);
yr = initialrect([3 3 4 4 3]);

him = findobj(gca, 'Type','image');
set(him,'HitTest','off');

hbox = line(xr,yr, 'Color',opt.color, 'LineStyle',':');

xresize = initialrect(1) + [0   w/2 w   w    w   w/2 0   0];
yresize = initialrect(3) + [0   0   0   h/2  h   h   h   h/2];

hresize = line(xresize,yresize, 'LineStyle','none', ...
    'Marker','s', 'MarkerFaceColor',opt.color, 'MarkerEdgeColor','none', ...
    'ButtonDownFcn',@resizerectangle);

hand = [hbox; hresize];

set(hbox, 'ButtonDownFcn',@(h,e) (dragrectangle(h,e, hand)));
set(hresize, 'ButtonDownFcn',@(h,e) (resizerectangle(h,e, hand)));

set(gcf, 'KeyPressFcn', @keypress, 'UserData',0);

waitfor(gcf,'UserData');

if (get(gcf,'UserData') == 1)
    xr = get(hbox,'XData');
    yr = get(hbox,'YData');
    
    finalrect = [xr(1) xr(2) yr(1) yr(3)];
else
    finalrect = [];
end;

delete(hand);

function dragrectangle(hobj, eventdata, hand)

pt0 = get(gca,'CurrentPoint');
set(gca, 'UserData',pt0);

for i = 1:length(hand),
    xd = get(hand(i),'XData');
    yd = get(hand(i),'YData');
    set(hand(i),'UserData',[xd(:)'; yd(:)']);
end;

set(gcf,'WindowButtonMotionFcn',@(h,e) (dragrectanglemove(h,e, hand)));
set(gcf,'WindowButtonUpFcn',@(h,e) (dragrectanglestop(h,e, hand)));


function dragrectanglemove(hobj, eventdata, hand)

pt1 = get(gca, 'CurrentPoint');
pt0 = get(gca, 'UserData');

dx = pt1(1,1) - pt0(1,1);
dy = pt1(1,2) - pt0(1,2);

for i = 1:length(hand)
    xy = get(hand(i),'UserData');
    set(hand(i),'XData',xy(1,:)+dx, 'YData',xy(2,:)+dy);
end;

function dragrectanglestop(hobj, eventdata, hand)

dragrectanglemove(hobj,eventdata, hand);
set(gcf,'WindowButtonMotionFcn',[], 'WindowButtonUpFcn',[]);







function resizerectangle(hobj, eventdata, hand)

pt0 = get(gca,'CurrentPoint');
xr = get(hand(2),'XData');
yr = get(hand(2),'YData');

setx = {[0 0 0 0 0], [0 0 0 0 0 0 0 0]};
sety = {[0 0 0 0 0], [0 0 0 0 0 0 0 0]};

[~,ind] = min(abs(pt0(1,1)-xr) + abs(pt0(1,2)-yr));
switch ind,
    case 1,
        setx = {[1 0 0 1 1], [1 0 0 0 0 0 1 1]};
        sety = {[1 1 0 0 1], [1 1 1 0 0 0 0 0]};

    case 2,
        sety = {[1 1 0 0 1], [1 1 1 0 0 0 0 0]};
        
    case 3,
        setx = {[0 1 1 0 0], [0 0 1 1 1 0 0 0]};
        sety = {[1 1 0 0 1], [1 1 1 0 0 0 0 0]};
        
    case 4,
        setx = {[0 1 1 0 0], [0 0 1 1 1 0 0 0]};
        
    case 5,
        setx = {[0 1 1 0 0], [0 0 1 1 1 0 0 0]};
        sety = {[0 0 1 1 0], [0 0 0 0 1 1 1 0]};
        
    case 6,
        sety = {[0 0 1 1 0], [0 0 0 0 1 1 1 0]};
        
    case 7,
        setx = {[1 0 0 1 1], [1 0 0 0 0 0 1 1]};
        sety = {[0 0 1 1 0], [0 0 0 0 1 1 1 0]};
        
    case 8,
        setx = {[1 0 0 1 1], [1 0 0 0 0 0 1 1]};
end;

set(gcf,'WindowButtonMotionFcn',@(h,e) (resizerectanglemove(h,e, hand,setx,sety)));
set(gcf,'WindowButtonUpFcn',@(h,e) (resizerectanglestop(h,e, hand,setx,sety)));


function resizerectanglemove(hobj, eventdata, hand,setx,sety)

pt1 = get(gca, 'CurrentPoint');

for i = 1:length(hand)
    x = get(hand(i),'XData');
    y = get(hand(i),'YData');
    x(setx{i} == 1) = pt1(1,1);
    y(sety{i} == 1) = pt1(1,2);
    
    set(hand(i),'XData',x, 'YData',y);
end;

function resizerectanglestop(hobj, eventdata, hand,setx,sety)

resizerectanglemove(hobj,eventdata, hand,setx,sety);
set(gcf,'WindowButtonMotionFcn',[], 'WindowButtonUpFcn',[]);




function keypress(hobj, eventdata)

c = get(hobj, 'CurrentCharacter');

switch (c),
    case {'q',char(10)}
        set(hobj,'UserData',1);
        
    case {'c',char(27)},
        set(hobj,'UserData',-1);
end;
