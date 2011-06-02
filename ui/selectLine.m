function [x,y,butt] = selectLine(param)

global HLINE1 HLINE2 LINEX LINEY KEY;

if (nargin > 0)
	eval(param);
	return;
end;

PT = 1;
LINEX = [];
LINEY = [];
KEY = [];

oldFcn = get(gcf, 'WindowButtonDownFcn');
oldDB = get(gcf, 'DoubleBuffer');
oldKey = get(gcf, 'KeyPressFcn');

HLINE1 = line('XData',LINEX, 'YData',LINEY, 'Visible','off',...
			 'Clipping','off','Color','k','LineStyle','-','UserData','');
HLINE2 = line('XData',LINEX, 'YData',LINEY, 'Visible','off',...
			 'Clipping','off','Color','w','LineStyle',':','UserData','');
			 
set(gcf, 'WindowButtonDownFcn', 'selectLine(''ButtonDown'');', ...
		 'KeyPressFcn', 'selectLine(''KeyPress'');',...
		 'DoubleBuffer', 'on');
		 
waitfor(HLINE1, 'UserData','done');
delete(HLINE1);
delete(HLINE2);

set(gcf, 'WindowButtonDownFcn', oldFcn, ...
		 'DoubleBuffer', oldDB,...
		 'KeyPressFcn', oldKey);

if (ischar(KEY))
	x = [];
	y = [];
	butt = KEY;
else
	x = LINEX;
	y = LINEY;
	butt = KEY;
end;

clear global HLINE1 HLINE2 LINEX LINEY;

%%%%%%%%%%%%%
%% KeyPress
function KeyPress

global HLINE1 HLINE2 KEY;

KEY = get(gcf, 'CurrentCharacter');
set(HLINE1, 'UserData','done');

%%%%%%%%%%%%%
%% ButtonDown
function ButtonDown

global HLINE1 HLINE2 LINEX LINEY PT KEY;

a = get(gca, 'CurrentPoint');
LINEX(1) = a(1,1);
LINEY(1) = a(1,2);

switch (get(gcf, 'SelectionType'))
case 'normal',
	KEY = 1;
case {'extend','alt'},
	KEY = 2;
case 'open',
	return;
end;

set(HLINE1, 'XData',LINEX, 'YData',LINEY, 'Visible','on');
set(HLINE2, 'XData',LINEX, 'YData',LINEY, 'Visible','on');
set(gcf, 'WindowButtonMotionFcn','selectLine(''MouseMove'');',...
		 'WindowButtonUpFcn', 'selectLine(''ButtonUp'');');

%%%%%%%%%%%%%
%% MouseMove
function MouseMove

global HLINE1 HLINE2 LINEX LINEY PT;

a = get(gca, 'CurrentPoint');
LINEX(2) = a(1,1);
LINEY(2) = a(1,2);

set(HLINE1, 'XData',LINEX, 'YData',LINEY, 'Visible','on');
set(HLINE2, 'XData',LINEX, 'YData',LINEY, 'Visible','on');

%%%%%%%%%%%%%
%% ButtonUp
function ButtonUp

global HLINE1 HLINE2 LINEX LINEY;

a = get(gca, 'CurrentPoint');
if ((abs(a(1,1)-LINEX(1)) > 1e-5) & (abs(a(1,2)-LINEY(1)) > 1e-5))	% different point
	LINEX(2) = a(1,1);
	LINEY(2) = a(1,2);
end;

set(HLINE1, 'XData',LINEX, 'YData',LINEY, 'Visible','on', 'UserData','done');
set(HLINE2, 'XData',LINEX, 'YData',LINEY, 'Visible','on');
set(gcf, 'WindowButtonMotionFcn','', 'WindowButtonUpFcn', '');
