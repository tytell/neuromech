function [x,y,r,butt] = selectLine(param)

global HCIRC1 HCIRC2 CIRCX CIRCY CTRX CTRY RADIUS KEY;
global THETA;

if (nargin > 0)
	eval(param);
	return;
end;

PT = 1;
CIRCX = [];
CIRCY = [];
CTRX = [];
CTRY = [];
RADIUS = 0;
KEY = [];

THETA = linspace(0,2*pi,20);

oldFcn = get(gcf, 'WindowButtonDownFcn');
oldDB = get(gcf, 'DoubleBuffer');
oldKey = get(gcf, 'KeyPressFcn');

HCIRC1 = line('XData',CIRCX, 'YData',CIRCY, 'Visible','off',...
			 'Clipping','off','Color','k','LineStyle','-','UserData','');
HCIRC2 = line('XData',CIRCX, 'YData',CIRCY, 'Visible','off',...
			 'Clipping','off','Color','w','LineStyle',':','UserData','');
			 
set(gcf, 'WindowButtonDownFcn', 'selectCircle(''ButtonDown'');', ...
		 'KeyPressFcn', 'selectCircle(''KeyPress'');',...
		 'DoubleBuffer', 'on');
		 
waitfor(HCIRC1, 'UserData','done');
delete(HCIRC1);
delete(HCIRC2);

set(gcf, 'WindowButtonDownFcn', oldFcn, ...
		 'DoubleBuffer', oldDB,...
		 'KeyPressFcn', oldKey);

if (ischar(KEY))
	x = [];
	y = [];
	butt = KEY;
else
	x = CTRX;
	y = CTRY;
	r = RADIUS;
	butt = KEY;
end;

clear global HCIRC1 HCIRC2 CIRCX CIRCY;

%%%%%%%%%%%%%
%% KeyPress
function KeyPress

global HCIRC1 HCIRC2 KEY;

KEY = get(gcf, 'CurrentCharacter');
set(HCIRC1, 'UserData','done');

%%%%%%%%%%%%%
%% ButtonDown
function ButtonDown

global CTRX CTRY PT KEY;

a = get(gca, 'CurrentPoint');
CTRX = a(1,1);
CTRY = a(1,2);

switch (get(gcf, 'SelectionType'))
case 'normal',
	KEY = 1;
case {'extend','alt'},
	KEY = 2;
case 'open',
	return;
end;

set(gcf, 'WindowButtonMotionFcn','selectCircle(''MouseMove'');',...
		 'WindowButtonUpFcn', 'selectCircle(''ButtonUp'');');

%%%%%%%%%%%%%
%% MouseMove
function MouseMove

global HCIRC1 HCIRC2 CIRCX CIRCY CTRX CTRY RADIUS THETA PT;

a = get(gca, 'CurrentPoint');
RADIUS = sqrt((a(1,1) - CTRX)^2 + (a(1,2) - CTRY)^2);

CIRCX = RADIUS * cos(THETA) + CTRX;
CIRCY = RADIUS * sin(THETA) + CTRY;

set(HCIRC1, 'XData',CIRCX, 'YData',CIRCY, 'Visible','on');
set(HCIRC2, 'XData',CIRCX, 'YData',CIRCY, 'Visible','on');

%%%%%%%%%%%%%
%% ButtonUp
function ButtonUp

global HCIRC1 HCIRC2 CTRX CTRY RADIUS;

a = get(gca, 'CurrentPoint');
if ((abs(a(1,1)-CTRX) < 1e-5) & (abs(a(1,2)-CTRY) < 1e-5))	% different point
	RADIUS = 0;
end;

set(HCIRC1, 'Visible','on', 'UserData','done');
set(HCIRC2, 'Visible','on');
set(gcf, 'WindowButtonMotionFcn','', 'WindowButtonUpFcn', '');
