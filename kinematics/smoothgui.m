function [sp,MSE,actMSE] = smoothgui(varargin)
% SMOOTHGUI Application M-file for smoothGui.fig
%    FIG = SMOOTHGUI launch smoothGui GUI.
%    SMOOTHGUI('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.5 18-Sep-2006 10:10:50

global SGt SGx SGy SGsp SGspx SGspy;
global SGxs SGys SGdxs SGdys SGddxs SGddys;
global SGxMSE SGyMSE;
global SGisJoint;

if (~isempty(varargin{1})),  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	
	set(handles.xMSESlider,'Min',0,'Max',5,'Value',0.5);
	set(handles.yMSESlider,'Min',0,'Max',5,'Value',0.5);
	set(handles.JointMSECheck,'Value',0);
	set(handles.ShowXCheck,'Value',1);
	set(handles.ShowYCheck,'Value',1);
	
	SGt = varargin{1};
	SGx = varargin{2};
	SGy = varargin{3};
	SGisJoint = 0;

	if (any(size(SGt) ~= size(SGx)) | any(size(SGx) ~= size(SGy))),
		error('x, y, and t must be the same size.');
	end;
	
	Smooth(handles);
	Draw(handles);
	
	guidata(fig, handles);
	
	waitfor(fig);
	
	if (SGisJoint),
		sp = SGsp;
		MSE = SGxMSE^2;
	else
		sp = [SGspx SGspy];
		MSE = [SGxMSE SGyMSE].^2;
	end;
	
	actxMSE = 1/sum(isfinite(SGx)) * nansum((SGx-SGxs).^2);
	actyMSE = 1/sum(isfinite(SGy)) * nansum((SGy-SGys).^2);
	actMSE = [actxMSE actyMSE];
	
elseif ischar(varargin{2}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

	try
		feval(varargin{2:end}); % FEVAL switchyard
	catch
		disp(lasterr);
	end;
end;

function Smooth(handles)

global SGt SGx SGy SGsp SGspx SGspy;
global SGxs SGys SGdxs SGdys SGddxs SGddys;
global SGxMSE SGyMSE;
global SGisJoint;

SGxMSE = get(handles.xMSESlider,'Value');
SGyMSE = get(handles.yMSESlider,'Value');

if (SGisJoint),
	k = find(isfinite(SGx) & isfinite(SGy));
	SGsp = spaps(SGt(k), [SGx(k); SGy(k)], SGxMSE^2, 3, ones(1,length(k))/length(k));
	
	xy = fnval(SGsp, SGt(k(1):k(end)));
	SGxs = repmat(NaN,size(SGx));
	SGxs(k(1):k(end)) = xy(1,:);
	SGys = repmat(NaN,size(SGy));
	SGys(k(1):k(end)) = xy(2,:);

	dxy = fnval(fnder(SGsp,1), SGt(k(1):k(end)));
	SGdxs = repmat(NaN,size(SGx));
	SGdxs(k(1):k(end)) = dxy(1,:);
	SGdys = repmat(NaN,size(SGy));
	SGdys(k(1):k(end)) = dxy(2,:);

	ddxy = fnval(fnder(SGsp,2), SGt(k(1):k(end)));
	SGddxs = repmat(NaN,size(SGx));
	SGddxs(k(1):k(end)) = ddxy(1,:);
	SGddys = repmat(NaN,size(SGy));
	SGddys(k(1):k(end)) = ddxy(2,:);
else
	k = find(isfinite(SGx) & isfinite(SGy));
	SGspx = spaps(SGt(k), SGx(k), SGxMSE^2, 3, ones(1,length(k))/length(k));
	SGspy = spaps(SGt(k), SGy(k), SGyMSE^2, 3, ones(1,length(k))/length(k));
	
	SGxs = repmat(NaN,size(SGx));
	SGxs(k(1):k(end)) = fnval(SGspx, SGt(k(1):k(end)));
	SGys = repmat(NaN,size(SGy));
	SGys(k(1):k(end)) = fnval(SGspy, SGt(k(1):k(end)));
	
	SGdxs = repmat(NaN,size(SGx));
	SGdxs(k(1):k(end)) = fnval(fnder(SGspx,1), SGt(k(1):k(end)));
	SGdys = repmat(NaN,size(SGy));
	SGdys(k(1):k(end)) = fnval(fnder(SGspy,1), SGt(k(1):k(end)));
	
	SGddxs = repmat(NaN,size(SGx));
	SGddxs(k(1):k(end)) = fnval(fnder(SGspx,2), SGt(k(1):k(end)));
	SGddys = repmat(NaN,size(SGy));
	SGddys(k(1):k(end)) = fnval(fnder(SGspy,2), SGt(k(1):k(end)));
end;

function Draw(handles)

global SGt SGx SGy SGsp SGspx SGspy;
global SGxs SGys SGdxs SGdys SGddxs SGddys;
global SGxMSE SGyMSE;
global SGisJoint;

axes(handles.PositionAxes);

x0 = nanmean(SGx);
y0 = nanmean(SGy);

leg = {};
if (get(handles.ShowXCheck,'Value')),
	plot(SGt, SGx - x0, 'r.', SGt, SGxs - x0, 'k-');
    leg = {'X','Smoothed X'};
	hold on;
end;
if (get(handles.ShowYCheck,'Value')),
	plot(SGt, SGy - y0, 'b.', SGt, SGys - y0, 'k-');
    leg = {leg{:}, 'Y','Smoothed Y'};
end;
hold off;
axis tight;

set(gca,'XTickLabel',{});
ylabel('Position (mm)');
legend(leg{:},'Location','Best');

axes(handles.VelAxes);

if (get(handles.ShowXCheck,'Value')),
	plot(SGt, SGdxs, 'r-');
	hold on;
end;
if (get(handles.ShowYCheck,'Value')),
	plot(SGt, SGdys, 'b-');
end;
hold off;
axis tight;

set(gca,'XTickLabel',{});
ylabel('Velocity (mm/s)');

axes(handles.AccelAxes);
if (get(handles.ShowXCheck,'Value')),
	plot(SGt, SGddxs, 'r-');
	hold on;
end;
if (get(handles.ShowYCheck,'Value')),
	plot(SGt, SGddys, 'b-');
end;
hold off;
axis tight;

ylabel('Acceleration (mm^2/s)');
xlabeL('Time (s)');

%| ABOUT CALLBACKS:
%| GUIDE automatically appends subfunction prototypes to this file, and 
%| sets objects' callback properties to call them through the FEVAL 
%| switchyard above. This comment describes that mechanism.
%|
%| Each callback subfunction declaration has the following form:
%| <SUBFUNCTION_NAME>(H, EVENTDATA, HANDLES, VARARGIN)
%|
%| The subfunction name is composed using the object's Tag and the 
%| callback type separated by '_', e.g. 'slider2_Callback',
%| 'figure1_CloseRequestFcn', 'axis1_ButtondownFcn'.
%|
%| H is the callback object's handle (obtained using GCBO).
%|
%| EVENTDATA is empty, but reserved for future use.
%|
%| HANDLES is a structure containing handles of components in GUI using
%| tags as fieldnames, e.g. handles.figure1, handles.slider2. This
%| structure is created at GUI startup using GUIHANDLES and stored in
%| the figure's application data using GUIDATA. A copy of the structure
%| is passed to each callback.  You can store additional information in
%| this structure at GUI startup, and you can change the structure
%| during callbacks.  Call guidata(h, handles) after changing your
%| copy to replace the stored original so that subsequent callbacks see
%| the updates. Type "help guihandles" and "help guidata" for more
%| information.
%|
%| VARARGIN contains any extra arguments you have passed to the
%| callback. Specify the extra arguments by editing the callback
%| property in the inspector. By default, GUIDE sets the property to:
%| <MFILENAME>('<SUBFUNCTION_NAME>', gcbo, [], guidata(gcbo))
%| Add any extra arguments after the last argument, before the final
%| closing parenthesis.



% --------------------------------------------------------------------
function xMSESlider_Callback(h, eventdata, handles, varargin)

set(handles.xMSEEdit, 'String', num2str(get(h,'Value')));
Smooth(handles);
Draw(handles);

guidata(handles.Figure, handles);

% --------------------------------------------------------------------
function xMSEEdit_Callback(h, eventdata, handles, varargin)

v = str2num(get(h,'String'));
if (~isempty(v)),
	set(handles.xMSESlider, 'Value', v);
	Smooth(handles);
	Draw(handles);
	
	guidata(handles.Figure, handles);
end;

% --------------------------------------------------------------------
function yMSESlider_Callback(h, eventdata, handles, varargin)

set(handles.yMSEEdit, 'String', num2str(get(h,'Value')));
Smooth(handles);
Draw(handles);

guidata(handles.Figure, handles);

% --------------------------------------------------------------------
function yMSEEdit_Callback(h, eventdata, handles, varargin)

v = str2num(get(h,'String'));
if (~isempty(v)),
	set(handles.yMSESlider, 'Value', v);
	Smooth(handles);
	Draw(handles);
	
	guidata(handles.Figure, handles);
end;

% --------------------------------------------------------------------
function JointMSECheck_Callback(h, eventdata, handles, varargin)

global SGisJoint;

if (get(h,'Value')),
	set(handles.yMSESlider,'Enable','off');
	set(handles.yMSEEdit,'Enable','off');
	
	SGisJoint = 1;
else
	set(handles.yMSESlider,'Enable','on');
	set(handles.yMSEEdit,'Enable','on');
	
	SGisJoint = 0;
end;

Smooth(handles);
Draw(handles);

guidata(handles.Figure, handles);

% --------------------------------------------------------------------
function ShowXCheck_Callback(h, eventdata, handles, varargin)

Draw(handles);


% --------------------------------------------------------------------
function ShowYCheck_Callback(h, eventdata, handles, varargin)

Draw(handles);


% --- Executes on button press in OKButton.
function OKButton_Callback(hObject, eventdata, handles)

close(handles.Figure);
