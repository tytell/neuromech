function thresh = thresholddlg(t,data, varargin)
% function thresh = thresholddlg(t,data, varargin)
% Shows a dialog to select a threshold on a 1D data set.
%
% Options:
%   'default', thresh - Default threshold value
%   'threshtype', ('abs','max', or 'min') - Type of threshold.  'abs'
%      selects values with magnitudes above the threshold, 'max' selects
%      raw value above the threshold, and 'min' selects raw value below the
%      threshold.
%   'name', str - Name of the dialog
%
% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2010, Eric Tytell <tytell at jhu dot edu>


opt.default = [];
opt.threshtype = 'abs';
opt.name = 'Click to set threshold';
opt.downsample = 1;

if ((ndims(data) > 2) || ~isvector(data))
    error('thresholddlg:VectorInputRequired','Data must be a vector');
end;

opt = parsevarargin(opt,varargin, 3, ...
    'multival',{'threshtype',{'abs','max','min'}});

dlg = dialog('Name',opt.name,...
    'WindowStyle','normal','NumberTitle','off', 'Units','pixels', ...
    'ToolBar','figure','Resize','on');
ax = axes('Position',[0.1 0.2 0.8 0.7], 'Parent',dlg);

hzm = zoom(dlg);
set(hzm, 'Motion','horizontal');

pos = get(dlg,'Position');

uicontrol('Style','pushbutton', 'String','OK', ...
    'Units','pixels','Position',[20 20 70 23],'Callback',@OkButtonClick, ...
    'Parent',dlg);
uicontrol('Style','pushbutton', 'String','Cancel', ...
    'Units','pixels','Position',[pos(3)-70-20 20 70 23], ...
    'Callback',@CancelButtonClick, 'Parent',dlg);

if (opt.downsample > 1)
    len = floor(length(data)/opt.downsample);
    data = reshape(data(1:len*opt.downsample),[opt.downsample len]);
    switch opt.threshtype,
        case 'abs',
            data = [min(data); max(data)];
        case 'max',
            data = max(data);
    end;
    data = data(:);
end;

if (numel(t) == 1)
    dt = t * opt.downsample;
    t = (0:length(data)-1)*dt;
elseif (opt.downsample > 1),
    d2 = ceil(opt.downsample/2);
    t = t(d2:opt.downsample:end);
end;

h = plot(ax,t,data,'k-', t,data,'g-', 'HitTest','off');
yl = prctile(data,[0.1 99.9]);
set(ax,'YLim',yl);

if (isempty(opt.default)),
    opt.default = prctile(data,85);
end;

switch lower(opt.threshtype),
    case 'abs',
        h(3) = line('XData',[t(1); t(end)],'YData',[opt.default; opt.default],...
            'Color','r', 'Parent',ax);
        h(4) = line('XData',[t(1); t(end)],'YData',-[opt.default; opt.default],...
            'Color','r', 'Parent',ax);
        
        isthresh = abs(data) <= opt.default;
        data1 = data;
        data1(isthresh) = NaN;
        set(h(2),'YData',data1);
        
        set(h,'HitTest','off');
        set(ax,'ButtonDownFcn',{@setThresh,h(1:4)});
        set(dlg,'KeyPressFcn',@keyPress, 'UserData','');
    case 'max',
        h(3) = line('XData',[t(1); t(end)],'YData',[opt.default; opt.default],...
            'Color','r', 'Parent',ax);
        
        isthresh = data <= opt.default;
        data1 = data;
        data1(isthresh) = NaN;
        set(h(2),'YData',data1);
        
        set(h,'HitTest','off');
        set(ax,'ButtonDownFcn',{@setThresh,h(1:3)});
        set(dlg,'KeyPressFcn',@keyPress, 'UserData','');
end;

waitfor(dlg,'UserData');

if (ishandle(dlg)),
    t = get(dlg,'UserData');
    switch t,
        case 'ok',
            thresh = get(h(3),'YData');
            thresh = thresh(1);
        case 'cancel',
            thresh = NaN;
    end;
    close(dlg);
else
    thresh = NaN;
end;

%------------------------------------------------------------------------
function setThresh(obj,eventdata, hLine)

pt = get(obj,'CurrentPoint');

thresh = pt(1,2);
data1 = get(hLine(1),'YData');
if (length(hLine) == 4)
    thresh = abs(thresh);
    set(hLine(3),'YData',[thresh; thresh]);
    set(hLine(4),'YData',-[thresh; thresh]);
    data1(abs(data1) <= thresh) = NaN;
elseif (length(hLine) == 3)
    set(hLine(3),'YData',[thresh; thresh]);
    data1(data1 <= thresh) = NaN;
end;

set(hLine(2),'YData',data1);


%------------------------------------------------------------------------
function keyPress(obj,eventdata)

c = get(obj, 'CurrentCharacter');
if (c == char(13)),                 % return
    set(obj,'UserData','ok');
end;

%------------------------------------------------------------------------
function OkButtonClick(obj,eventdata)

fig = get(obj,'Parent');
set(fig,'UserData','ok');


%------------------------------------------------------------------------
function CancelButtonClick(obj,eventdata)

fig = get(obj,'Parent');
set(fig,'UserData','cancel');

