function data = analyzePIVvortices(data)

% first pull up the dialog to get the frames to operate on
% and the threshold parameters and so forth
dlg = openfig(mfilename, 'new');
dlgdata = guihandles(dlg);

sz = size(data.PIV.u);
if (length(sz) == 2),
    sz(3) = 1;
end;
dlgdata.Figure = dlg;
dlgdata = initializeGUI(dlgdata, sz);

guidata(dlg,dlgdata);

uiwait(dlg);

% if everything went well, the figure is still around
if (ishandle(dlg)),
    dlgdata = guidata(dlg);
    delete(dlg);

    % points on a unit circle for plotting purposes
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

    if (dlgdata.Frames(1) < 1),
        dlgdata.Frames(1) = 1;
    end;
    if (dlgdata.Frames(2) > data.nFrames)
        dlgdata.Frames(2) = data.nFrames;
    end;

    % convert the velocity length unit into the position unit
    % since identifyVortices returns positions (mostly), we don't want
    % to convert position units, because then identifyVortices will return
    % the positions in the converted units, which is confusing
    [velscalefac,unitserror] = feval(data.generalFcns.apConvertUnits, ...
                                     data.Units.vel{1}, data.Units.pos);

    % actually do the vortex identification
    % unfortunately, there's no good way to do a timed wait bar
    % so we'll just display a "please wait" box
    k = dlgdata.Frames(1):dlgdata.Frames(2);
    vxr = identifyVortices(data.PIV.x,data.PIV.y,...
                           data.PIV.u(:,:,k)*velscalefac,...
                           data.PIV.v(:,:,k)*velscalefac,...
                           'threshold',dlgdata.Thresh,...
                           'minframes',dlgdata.MinFrames,...
                           'nocirc');

    % translate the vxr structure to the region structure
    norig = size(data.Regions,1);
    for i = 1:length(vxr),
        cur = norig+i;
        data.Regions(cur,1) = struct('handle',-1, 'type','circle', ...
                                     'label',-1, 'x',[],'y',[],'status',0);
        data.Regions(cur,2:data.nFrames) = ...
            repmat(data.Regions(cur,1),[1 data.nFrames-1]);

        for j = 1:length(vxr(i).frames),
            data.Regions(cur,vxr(i).frames(j)).x = vxr(i).ctrx(j) + ...
                1.5*vxr(i).rgeom(j)*x0;
            data.Regions(cur,vxr(i).frames(j)).y = vxr(i).ctry(j) + ...
                1.5*vxr(i).rgeom(j)*y0;
            data.Regions(cur,vxr(i).frames(j)).status = 2;
        end;
    end;

    % all of the regions need to be plotted in the current frame, even
    % if they don't have points
    figure(data.Figure);
    rgns = data.Regions(:,data.curFrame);
    lastgood = NaN;
    for i = 1:length(vxr),
        cur = norig+i;

        if (isempty(rgns(cur).x)),
            viz = 'off';
            xl = 0;
            yl = 0;
        else
            viz = 'on';
            xl = rgns(cur).x(end);
            yl = rgns(cur).y(end);
            lastgood = cur;
        end;

        rgns(cur).handle = ...
            line('XData',rgns(cur).x, 'YData',rgns(cur).y, ...
                 'EraseMode','normal','Color','k',...
                 'ButtonDownFcn', {data.rgnFcns.apClickRgn, data.Panel}, ...
                 'Tag','Vx','Visible',viz);
        data.rgnNames{cur} = sprintf('Vx %d',data.rgnNameNum+i);
        rgns(cur).label = ...
            text(xl,yl, data.rgnNames{cur},...
                 'HorizontalAlignment','left',...
                 'Color','k','HitTest','off','Visible',viz);
    end;
    data.Regions(:,data.curFrame) = rgns;
    data.rgnNameNum = data.rgnNameNum + length(vxr);

    data = feval(data.calcGuiFcns.apAddCalc, data);
    data = feval(data.rgnFcns.apSetCurRgn, data, lastgood);
end;


% -----------------------------------------------
function GoButton_Callback(obj,eventdata)

data = guidata(obj);

if (get(data.AllFramesCheck,'Value')),
    data.Frames = [1 Inf];
else
    data.Frames(1) = str2num(get(data.FrameStartEdit,'String'));
    data.Frames(2) = str2num(get(data.FrameEndEdit,'String'));
end;

data.Thresh = str2num(get(data.DCEVThreshEdit,'String'));
data.MinFrames = str2num(get(data.MinDurationEdit,'String'));

guidata(obj,data);
uiresume;

% -----------------------------------------------
function CancelButton_Callback(obj,eventdata)

dlg = get(obj,'Parent');
delete(dlg);

% -----------------------------------------------
function data = initializeGUI(data, sz)

set(data.FrameStartEdit, 'String','1');
set(data.FrameEndEdit, 'String',num2str(sz(3)));
set(data.DCEVThreshEdit, 'String','300');
set(data.MinDurationEdit, 'String','1');
set(data.GoButton, 'Callback', @GoButton_Callback);
set(data.CancelButton, 'Callback', @CancelButton_Callback);

                  
