function varargout = analyzepiv(varargin)
% analyzepiv - Prompts to load data from DAT files or a Matlab file
%
% analyzepiv(filename) - Loads a Matlab file
%
% analyzepiv(x,y,t,u,v,[w],optionalparams) - Loads directly from the
%   variables passed in.  w is an optional "through-plane" z component of
%   velocity from stereo PIV.

constDEFAULT_UNITS = {'mm', 'm/s'};

if (nargin == 0),
    [piv,names,units] = analyzePIVload;
    if (isempty(piv)),
        return;
    end;
elseif ((nargin == 1) && ischar(varargin{1})),
    [piv,names,units] = analyzePIVload(varargin{1});
    if (isempty(piv)),
        return;
    end;
elseif (nargin >= 5),  % x,y,t,u,v, possibly w, then extras.
    param = varargin;

   %check for possible w velocity component
    if ((nargin >= 6) && (ndims(param{6}) == ndims(param{4})) && ...
        all(size(param{6}) == size(param{4}))),
        basenames = {'x','y','t','u','v','w'};
        paramInd = 6;
    else
        basenames = {'x','y','t','u','v'};
        paramInd = 5;
    end;

   %get any optional parameters, if they exist
    if (nargin > paramInd),
       %optional parameters have text names preceeding them (like
       %analyzePIV(..., 'vorticity',vort, ...)), so look for char
       %parameters.  We'll use the names to construct the PIV structure
        isch = find(cellfun('isclass',param,'char'));
        if ((max(isch)+1 > nargin) || ...
            any(~cellfun('isclass',param(isch+1),'double'))),
            error('All optional parameters must have numeric values.');
        end;

        optnames = param(isch);
        param = param(setdiff(1:nargin,isch));
    else
        optnames = {};
    end;

    piv = cell2struct(param,{basenames{:},optnames{:}},2);
    if (isfield(piv,'w')),
        piv.isStereo = true;
    else
        piv.w = [];
        piv.isStereo = false;
    end;
    piv.frames = 1:size(piv.u,3);
    piv.fileType = 'base';

    units.pos = '';
    units.time = '';
    units.vel = {'',''};
end;

if (any(isnan(piv.x(:)) | isnan(piv.y(:)))),
    % check for plaid-ness
    dx = diff(piv.x,[],2);
    dy = diff(piv.y,[],1);

    sz = size(piv.x);
    if (ndims(piv.x) == 2),
        sz(3) = 1;
    end;

    x = first(piv.x,isfinite(piv.x),1);
    dx0 = diff(x);
    y = first(piv.y,isfinite(piv.y),2);
    dy0 = diff(y);

    if (any(flatten(isfinite(dx) && ...
                    (dx ~= repmat(dx0,[sz(1) 1 sz(3)])))) || ...
        any(flatten(isfinite(dy) && ...
                    (dy ~= repmat(dy0,[1 sz(2:3)]))))),
        error('x and y must be plaid.');
    end;

    piv.x = repmat(x,[sz(1) 1 1]);
    piv.y = repmat(y,[1 sz(2) 1]);
end;

% set up the units
if (exist('units','var')),
    if (isfield(units,'x') && ~isempty(units.x)),
        u = {units.x, [units.vel{1} '/' units.vel{2}]};
    elseif (isfield(units,'pos') && ~isempty(units.pos)),
        u = {units.pos, [units.vel{1} '/' units.vel{2}]};
    else
        u = constDEFAULT_UNITS;  % defined at beginning of function
    end;
else
    u = constDEFAULT_UNITS;  % defined at beginning of function
end;

u = inputdlg({'Position units','Velocity units'},'Choose units',1,u);

units.pos = u{1};
tok = regexp(u{2},'(.+)/(.+)','tokens','once');
if (~isempty(tok)),
    units.vel{1} = tok{1};
    units.vel{2} = tok{2};
else
    units.vel = u{2};
end;

piv.smoothval = 0;

fig = openfig(mfilename, 'new');
data = guihandles(fig);

data.Panel = fig;
data.Figure = figure('WindowStyle','normal');
set(data.Figure,'UserData',fig);
data.PIV = piv;
data.SaveDataFile = '';
data.SaveDataType = '';
data.Units = units;
data.Background = '';
data.hBackground = -1;
data.subtractVector = [0 0 0];
data.Calc.doFcns = [];
data.Calc.nRows = 0;
data.Calc.TopRow = 1;
data.axisOrient = 'xy';

data.nFrames = size(data.PIV.u,3);
data.curFrame = 1;
data.frameSkip = 1;

data.Regions = struct('handle',{},'type',{},'label',{},...
                      'x',{},'y',{},'status',{});
data.rgnNames = {};
data.rgnSzHandle = -1;
data.rgnRotHandle = -1;

data.generalFcns.apConvertUnits = @apConvertUnits;
data.generalFcns.apSetFrame = @apSetFrame;

if (data.PIV.smoothval ~= 0)
    d = 0.5*((piv.x(1,2,1) - piv.x(1,1,1)) + (piv.y(2,1,1) - piv.y(1,1,1)));
    [piv.us,piv.vs] = agw(piv.x(:,:,1),piv.y(:,:,1), piv.u,piv.v, ...
        data.PIV.smoothval*d);
else
    piv.us = piv.u;
    piv.vs = piv.v;
end
%try to show 50x50 vectors
shw = 50/max([size(piv.u,1) size(piv.u,2)]);
if (shw > 1),
    shw = 1;
end;
shw = diground(shw^2,0.05);               % round to the nearest 5%
if (shw == 0)
    shw = 0.05;
end
[data.hVectors,data.vectorOpts] = quiverc(piv.x(:,:,1),piv.y(:,:,1),...
                        piv.us(:,:,1),piv.vs(:,:,1),'show',shw);
axis('equal','tight',data.axisOrient);
data.Axes = gca;
if (~isempty(data.Units.pos)),
    h = xlabel(strcat('x (',data.Units.pos,')'));
    ylabel(strcat('y (',data.Units.pos,')'));
else
    h = xlabel('x');
    ylabel('y');
end;
data.hUVText = text('Units','normalized','Position',[1 -0.1],...
                    'String','vec = (0.000, 0.000) [m/s]',...
                    'HorizontalAlignment','right');
data.OldMotionFcn = {@apVecMouseMove,data.Panel};
set(data.Figure,'WindowButtonMotionFcn',{@apVecMouseMove,data.Panel});

ax = [min(piv.x(:)) max(piv.x(:)) min(piv.y(:)) max(piv.y(:))];
rnguv = prctile(sqrt(piv.us(:).^2 + piv.vs(:).^2),95) * ...
        data.vectorOpts.Scale;
ax = ax + [-rnguv rnguv -rnguv rnguv];  % ensure vectors fit in the axes
axis(ax);
data.axisRange = ax;

data = apSetup(data);
guidata(fig,data);

uiwait(data.Panel);
if (ishandle(data.Panel)),
    data = guidata(data.Panel);
    delete(data.Panel);
    if (nargout == 1),
        varargout = {data.Regions};
    end;
else
    if (nargout == 1),
        varargout = {[]};
    end;
end;
if (ishandle(data.Figure)),
    set(data.Figure,'CloseRequestFcn','closereq', ...
                    'WindowButtonMotionFcn','');
    h = allchild(data.Axes);
    set(h,'ButtonDownFcn',[]);               % clear all callbacks
    if (ishandle(data.rgnSzHandle)),
        delete(data.rgnSzHandle);
    end;
    if (ishandle(data.rgnRotHandle)),
        delete(data.rgnRotHandle);
    end;
end;

% -------------------------------------------------
function apSave(obj,eventdata)

data = guidata(obj);

butt = questdlg('What information should be saved?','Save what?',...
                'Processed data','Regions','Regions & data',...
                'Processed data');

filt = {'*.csv','Excel compatible file (*.csv)';...
        '*.mat','Matlab data file (*.mat)'};

isdata = 0;
isrgns = 0;
switch butt,
 case 'Processed data',
  availfilt = 1:2;
  isdata = 1;
 case 'Regions',
  availfilt = 2;
  isrgns = 1;
 case 'Regions & data',
  availfilt = 2;
  isrgns = 1;
  isdata = 1;
end;

[file,path,filtind] = uiputfile(filt(availfilt,:),...
                        'Save...');
if (filtind == 0),
    % cancel button pressed
    return;
end;

file = fullfile(path,file);

if (isrgns),
    Regions = data.Regions;

    % check if the file exists - if it does, append the Regions struct
    if (exist(file,'file')),
        app = {'-append'};
    else
        app = {};
    end;
    save(file,'Regions',app{:});
end;

if (isdata),
    data.SaveDataFile = file;
    if ((length(availfilt) == 2) && (filtind == 1)),
        data.SaveDataType = 'csv';
    else
        data.SaveDataType = 'matlab';
    end;

    [q,data] = feval(data.calcGuiFcns.apSaveData,data);
end;

guidata(obj,data);


% -------------------------------------------------
function apSetFrame(obj,eventdata, type,frame)

data = guidata(obj);

prev = data.curFrame;

switch type,
 case 'slide',
  fr1 = round(get(obj,'Value'));

 case 'edit',
  fr1 = str2num(get(obj,'String'));
  
 case 'set',
  fr1 = frame;
end;

if (~isempty(fr1) && (fr1 >= 1) && (fr1 <= data.nFrames)),
    data.curFrame = round(fr1);
else
    beep;
end;

data = feval(data.rgnFcns.apSetRgnFrame,data,data.curFrame,prev);

set(data.FrameEdit,'String',num2str(data.curFrame));
set(data.FrameSlider,'Value',data.curFrame);

data = analyzePIVupdate(data, {'vectors','regions','background'},'recalc');
guidata(obj,data);

% -------------------------------------------------
function apSetFrameSkip(obj,eventdata)

data = guidata(obj);

skip = str2num(get(obj,'String'));
if (~isempty(skip) && (mod(skip,1) == 0) && ...
    (skip > 0) && (skip < data.nFrames)),
    data.frameSkip = round(skip);
    set(data.FrameSlider,'SliderStep',[skip/(data.nFrames-1) 0.1]);
end;

set(obj,'String',num2str(data.frameSkip));

guidata(obj,data);



% -------------------------------------------------
function apScale(obj,eventdata, type)

data = guidata(obj);

switch type,
 case 'slide',
  sg = get(obj,'Value');
  set(obj,'Value',0);
  if (sg > 0),                          % make longer
      sc = data.vectorOpts.Scale*1.25;
  else                                  % make shorter
      sc = data.vectorOpts.Scale/1.25;
  end;
 case 'set',
  sc1 = str2num(get(obj,'String'));
  if (~isempty(sc1) && (sc1 >= 0)),
      sc = sc1;
  end;
 case 'auto',
  sc = [];
end;

data.vectorOpts.AbsScale = sc;
data = analyzePIVupdate(data,{'vectors','regions'});

set(data.ScaleEdit,'String',num2str(data.vectorOpts.Scale));

guidata(obj,data);

% -------------------------------------------------
function apTruncate(obj,eventdata)

data = guidata(obj);

trunc = str2num(get(obj,'String'));
if (~isempty(trunc) && (trunc > 0) && (trunc <= 100)),
    data.vectorOpts.ScaleRange(2) = trunc/100;
    data = analyzePIVupdate(data,{'vectors','regions'},'rescale');
end;

set(obj,'String',num2str(round(data.vectorOpts.ScaleRange(2)*100)));

guidata(obj,data);

% -------------------------------------------------
function apSetShow(obj,eventdata)

data = guidata(obj);

shw = str2num(get(obj,'String'));
if (~isempty(shw) && (shw > 0) && (shw <= 100)),
    data.vectorOpts.Show = shw/100;
    data = analyzePIVupdate(data,{'vectors','regions'},'rescale');
end;

set(obj,'String',num2str(round(data.vectorOpts.Show*100)));

guidata(obj,data);

% -------------------------------------------------
function apSetHeadSize(obj,eventdata)

data = guidata(obj);

hsz = get(obj,'Value');
data.vectorOpts.HeadSize = hsz;
data = analyzePIVupdate(data,{'vectors','regions'});

guidata(obj,data);

% -------------------------------------------------
function apKeepHeads(obj,eventdata)

data = guidata(obj);

keep = str2num(get(obj,'String'));
if (~isempty(keep) && (keep > 0) && (keep <= 100)),
    data.vectorOpts.HeadRange(1) = keep/100;
    data.vectorOpts.AbsHeadRange = [];
    data = analyzePIVupdate(data,{'vectors','regions'});
end;

set(obj,'String',num2str(data.vectorOpts.HeadRange(1)*100));

guidata(obj,data);

% -------------------------------------------------
function apArrowColor(obj, eventdata)

data = guidata(obj);

strs = get(obj,'String');
ind = get(obj,'Value');

switch strs{ind},
 case 'Custom',
   % do nothing for the moment
 otherwise
  data = apSetArrowColor(data, strs(ind));
end;

guidata(obj,data);


% -------------------------------------------------
function [col,colstr,colstrnum] = apGetArrowColor(data, strs)

RGB = [0 0 1; 0 1 0; 1 0 0; 0 1 1; 1 0 1; 1 1 0; 1 1 1; 0 0 0];
Names = {'Blue'; 'Green'; 'Red'; 'Cyan'; 'Magenta'; 'Yellow'; ...
         'White'; 'Black'};

col = [];
colstr = '';
colstrnum = [];

switch (get(data.hVectors(1),'Type')),
 case 'line',
  col = get(data.hVectors(1),'Color');
 case 'patch',
  col = get(data.hVectors(1),'EdgeColor');
  if (~isnumeric(col)),
      col = get(data.hVectors(1),'FaceColor');
  end;
end;
if (isnumeric(col)),
    ind = find(all(RGB == repmat(col,[size(RGB,1) 1]),2));
    
    if (~isempty(ind)),
        colstr = Names{ind};
        ind2 = strmatch(lower(Names{ind}),lower(strs));
        if (~isempty(ind2)),
            colstrnum = ind2;
        end;
    end;
end;

% -------------------------------------------------
function data = apSetArrowColor(data, col)

setquivercol(data.hVectors,col);
data.vectorOpts.Color = col;

% -------------------------------------------------
function apBackgroundColor(obj, eventdata)

data = guidata(obj);

strs = get(obj,'String');
ind = get(obj,'Value');

switch strs{ind},
 case 'White',
  set(data.Axes,'Color','w');
  data.Background = '';
 case 'Black',
  set(data.Axes,'Color','k');
  data.Background = '';
 case 'Vorticity',
  data = apCalcDeriv(data,'vorticity');
  data.Background = 'vorticity';
  data.bgSymmetric = 1;

  data = analyzePIVupdate(data,'background');
 case 'Shear',
  data = apCalcDeriv(data,'shear');
  data.Background = 'shear';
  data.bgSymmetric = 1;

  data = analyzePIVupdate(data,'background');
 case 'DCEV',
  data = apCalcDeriv(data,'dcev');
  data.Background = 'dcev';
  data.bgSymmetric = 0;

  data = analyzePIVupdate(data,'background');
 case 'u',  % color by u component of velocity
  data.Background = 'u';
  data.bgSymmetric = 1;

  data = analyzePIVupdate(data,'background');
 case 'v',  % color by v component of velocity
  data.Background = 'v';
  data.bgSymmetric = 1;

  data = analyzePIVupdate(data,'background');
 case 'w',  % color by w component of velocity
  data.Background = 'w';
  data.bgSymmetric = 1;

  data = analyzePIVupdate(data,'background');
  
end;

guidata(obj,data);


% -------------------------------------------------
function data = apDoSmooth(data)

if (data.PIV.smoothval > 0)
    posscalefac = apConvertUnits(data.Units.pos, data.Units.vel{1});
    d = 0.5*((data.PIV.x(1,2,1) - data.PIV.x(1,1,1)) + (data.PIV.y(2,1,1) - data.PIV.y(1,1,1)));
    d = d*posscalefac;
    [data.PIV.us, data.PIV.vs, DU] = agw(posscalefac*data.PIV.x,posscalefac*data.PIV.y, ...
        data.PIV.u,data.PIV.v, ...
        posscalefac*data.PIV.x,posscalefac*data.PIV.y, data.PIV.smoothval*d);
else
    data.PIV.us = data.PIV.u;
    data.PIV.vs = data.PIV.v;
    DU = [];
end
if (ismember(data.Background,{'vorticity','shear','dcev'}))
    data.PIV.(data.Background) = [];
    data = apCalcDeriv(data,data.Background, DU);
end
data = analyzePIVupdate(data,{'vectors','background'});

% -------------------------------------------------
function apSmoothToggle(obj,eventdata)

data = guidata(obj);

issmooth = get(obj,'Value') > 0;

if (issmooth)
    set(data.SmoothEdit,'Enable','on');
    smoothval1 = str2double(get(data.SmoothEdit, 'String'));
    if (~isempty(smoothval1) && (smoothval1 > 0))
        data.PIV.smoothval = smoothval1;
    end
else
    set(data.SmoothEdit,'Enable','off');
    data.PIV.smoothval = 0;
end;
data = apDoSmooth(data);

guidata(obj,data);

% -------------------------------------------------
function apSmoothEdit(obj,eventdata)

data = guidata(obj);

smoothval1 = str2double(get(obj,'String'));
if (~isempty(smoothval1) && (smoothval1 > 0))
    data.PIV.smoothval = smoothval1;
    data = apDoSmooth(data);
end

guidata(obj,data);


% -------------------------------------------------
function data = apCalcDeriv(data, what, DU)

if (~isfield(data.PIV,what) || isempty(data.PIV.(what))),
    % need to rescale position so that it has the same units as
    % the length unit in velocity - otherwise the derivative ends up
    % weird
    posscalefac = apConvertUnits(data.Units.pos, data.Units.vel{1});

    if ((nargin == 2) || isempty(DU))
        if isfield(data.PIV,'us')
            DU = velderiv(posscalefac*data.PIV.x,posscalefac*data.PIV.y, ...
                data.PIV.us,data.PIV.vs);
        else
            DU = velderiv(posscalefac*data.PIV.x,posscalefac*data.PIV.y, ...
                data.PIV.u,data.PIV.v);
        end
    end
    
    switch what,
     case 'vorticity',
      data.PIV.vorticity = cat(3,DU.dvdx) - cat(3,DU.dudy);
     case 'shear',
      data.PIV.shear = cat(3,DU.dvdx) + cat(3,DU.dudy);
     case 'dcev',
      data.PIV.dcev = -((cat(3,DU.dudx) + cat(3,DU.dvdy)).^2 - ...
          4*(cat(3,DU.dudx) .* cat(3,DU.dvdy) - ...
             cat(3,DU.dudy) .* cat(3,DU.dvdx)));
      data.PIV.dcev(data.PIV.dcev < 0) = 0;
    end;
end;

% -------------------------------------------------
function [fac,converr] = apConvertUnits(from,to)

knownunits = {'m','cm','mm','um','micron','ft','in'};
conv2meters = [1 1e-2 1e-3 1e-6 1e-6 3.0480e-1 2.54e-2];

% if the units are the same, then it doesn't matter if we
% know them or not - the conversion is 1
if (strcmp(from,to)),
    fac = 1;
    converr = 0;
else
    % otherwise, try to sort out the conversion
    fromind = strmatch(lower(from),knownunits,'exact');
    toind = strmatch(lower(to),knownunits,'exact');

    if ((numel(fromind) == 1) && (numel(toind) == 1)),
        fac = conv2meters(fromind)/conv2meters(toind);
        converr = 0;
    else
        fac = 1;
        converr = 1;
    end;
end;

% -------------------------------------------------
function apSubtract(obj, eventdata)

data = guidata(obj);
bStereo = data.PIV.isStereo;

sub = get(obj,'Value');

switch sub,
 case 1,
  val = [0 0 0];
 case 2,  % subtract average flow
  val = nanmean([data.PIV.us(:) data.PIV.vs(:)]);
  if bStereo, val = [val, nanmean(data.PIV.ws(:))]; end
 case 3,  % substract average flow
  val = nanmedian([data.PIV.us(:) data.PIV.vs(:)]);
  if bStereo, val = [val, nanmean(data.PIV.w(:))]; end
 case 4, % subtract average u
  val = [nanmean(data.PIV.us(:)) 0 0];
 case 5, % subtract average v
  val = [0 nanmean(data.PIV.vs(:)) 0];
 case 6, % subtract average w
  if ~isempty(data.PIV.w), 
      val = [0 0 nanmean(data.PIV.ws(:))];
  else
      val = [0 0 0]; 
  end
 otherwise,
  warndlg('This feature is not supported yet.');
  val = [0 0 0];
end;

data.subtractVector = val;
data = analyzePIVupdate(data,'vectors');

guidata(obj,data);


% -------------------------------------------------
function apVecMouseMove(obj,eventdata,panel)

data = guidata(panel);

pt = get(data.Axes,'CurrentPoint');

if (size(data.PIV.x,3) > 1),
    px = data.PIV.x(1,:,data.curFrame);
    py = data.PIV.y(:,1,data.curFrame);
else
    px = data.PIV.x(1,:);
    py = data.PIV.y(:,1);
end;
[~,i] = min(abs(pt(1,2)-py));
[~,j] = min(abs(pt(1,1)-px));

if (isfield(data.PIV,'us'))
    u = data.PIV.us(i,j,data.curFrame) - data.subtractVector(1);
    v = data.PIV.vs(i,j,data.curFrame) - data.subtractVector(2);
else
    u = data.PIV.u(i,j,data.curFrame) - data.subtractVector(1);
    v = data.PIV.v(i,j,data.curFrame) - data.subtractVector(2);
end
if ~data.PIV.isStereo
    w = NaN;
elseif isfield(data.PIV,'ws')
    w = data.PIV.ws(i,j,data.curFrame) - data.subtractVector(3);
else
    w = data.PIV.w(i,j,data.curFrame) - data.subtractVector(3);
end

if (~isempty(data.Units.vel)),
    txt = sprintf('vec = (%8.3g,%.3g,%.3g) [%s/%s]',u,v,w, ...
                  data.Units.vel{1},data.Units.vel{2});
else
    txt = sprintf('vec = (%8.3g,%.3g,%.3g)',u,v,w);
end;

set(data.hUVText,'String',txt);

% note: no guidata(obj,data) for speed purposes


% -------------------------------------------------
function apCloseFigure(obj, eventdata, panel)

data = guidata(panel);
switch questdlg('Save data?','Quit','Yes'),
 case 'Yes',
  doclose = feval(data.calcGuiFcns.apSaveData,data);
 case 'No',
  doclose = 1;
 case 'Cancel',
  doclose = 0;
end;

if (doclose),
    delete(obj);
    if (ishandle(panel)),
        uiresume(panel);
    end;
end;

% -------------------------------------------------
function apKeyPress(obj, eventdata, panel)

data = guidata(panel);

c = get(obj,'CurrentCharacter');

switch lower(c),
 case {char(127), 'd'},
  action = questdlg('Delete region completely or stop after this frame?', ...
                    'Delete or stop?', 'Delete','Stop','Cancel', ...
                    'Cancel');
  switch action,
   case 'Delete',
    data = feval(data.rgnFcns.apDeleteRgn,data,data.curRgn);
   case 'Stop',
    data = feval(data.rgnFcns.apStopRgn,data,data.curRgn);
  end;
  
case {char(13),char(29)},                   % enter or right arrow
    apSetFrame(panel,eventdata, 'set',data.curFrame+data.frameSkip);
    data = guidata(panel);      % apSetFrame will have changed the data structure, so we update it here
    
case char(28),                              % left arrow
    apSetFrame(panel,eventdata, 'set',data.curFrame-data.frameSkip);
    data = guidata(panel);
    
end;

guidata(panel,data);

% -------------------------------------------------
function data = apSetup(data)

set(data.ScaleSlider,'Callback',{@apScale,'slide'});
set(data.ScaleEdit,'Callback',{@apScale,'set'},...
                  'String',num2str(data.vectorOpts.Scale));
set(data.ScaleAutoButton,'Callback',{@apScale,'auto'});
set(data.TruncateEdit,'Callback',@apTruncate, ...
         'String',num2str(round(data.vectorOpts.ScaleRange(2)*100)));
set(data.ShowButton,'Callback',@apSetShow, ...
                  'String',num2str(round(data.vectorOpts.Show*100)));
set(data.HeadSizeSlider,'Callback',@apSetHeadSize, ...
                  'Value',data.vectorOpts.HeadSize);
set(data.KeepHeadsEdit,'Callback',@apKeepHeads, ...
         'String',num2str(round(data.vectorOpts.HeadRange(1)*100)));

colstrs = get(data.ColorDrop,'String');
[col,colstr,colstrnum] = apGetArrowColor(data,colstrs);
set(data.ColorDrop,'Callback',@apArrowColor, ...
                  'Value',colstrnum);

set(data.BackgroundDrop,'Callback',@apBackgroundColor);
set(data.SubtractDrop,'Callback',@apSubtract);

set(data.SmoothCheck,'Callback',@apSmoothToggle);
set(data.SmoothEdit,'Callback',@apSmoothEdit);

[data,rgn] = analyzePIVrgn(data);
data.rgnFcns = rgn;

set(data.SaveButton,'Callback',@apSave);
set(data.SaveDataMenu,'Callback',@apSave);

[data,calcgui,calc] = analyzePIVcalc(data);
data.calcGuiFcns = calcgui;
data.calcFcns = calc;

if (data.nFrames > 1),
    set(data.FrameSlider,'Callback',{@apSetFrame,'slide'}, ...
                      'Min',1,'Max',data.nFrames,'Value',1, ...
                      'SliderStep',[1/(data.nFrames-1) 0.1]);
    set(data.FrameEdit,'Callback',{@apSetFrame,'edit'}, ...
                      'String','1');
    set(data.FrameSkipEdit,'Callback',@apSetFrameSkip, ...
                      'String','1');
else
    set(data.FrameSlider,'Enable','off');
    set(data.FrameEdit,'Enable','off');
    set(data.FrameSkipEdit,'Enable','off');
end;

set(data.Figure, 'KeyPressFcn', {@apKeyPress,data.Panel});
set(data.Figure, 'CloseRequestFcn', {@apCloseFigure,data.Panel});
set(data.Panel, 'CloseRequestFcn', {@apCloseFigure,data.Panel});
set(data.QuitMenu, 'Callback',{@apCloseFigure,data.Panel});


