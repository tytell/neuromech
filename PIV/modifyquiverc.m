function modifyquiverc(quiveraxes)

if (nargin == 0),
    quiveraxes = gca;
end;

hq = findobj(quiveraxes, 'Tag','quiverc');
if (isempty(hq)),
    error('No quiverc vectors in the current axes.');
elseif (length(hq) > 2),
    error('Cannot deal with multiple quiverc plots right now.');
end;
opts = get(hq(1),'UserData');

fig = openfig(mfilename, 'new');

data = guihandles(fig);
data.hQuiverc = hq;
data.Options = opts;

mqSetUICallbacks(data);
mqSetUIStart(data);

guidata(fig,data);


% -------------------------------------------------
function data = mqUpdateVectors(data)

delete(data.hQuiverc);
rep = get(data.Options.Axes,'NextPlot');
set(data.Options.Axes,'NextPlot','add');
[data.hQuiverc,data.Options] = quiverc(data.Options);
set(data.Options.Axes,'NextPlot',rep);

% -------------------------------------------------
function mqSetColor(obj,eventdata, type)

data = guidata(obj);

switch type,
 case 'name',
  strs = get(obj,'String');
  num = get(obj,'Value');
  col = strs{num};
 case 'RGB',
  str = get(obj,'String');
  if (str(1) == '['),
      str = str(2:end);
  end;
  col = sscanf(str,'%f')';
  if ((length(col) ~= 3) | any(col < 0) | any(col > 1)),
      beep;
      set(obj,'String','');
  end;
end;

ec = get(data.hQuiverc,'EdgeColor');
i = find(~cellfun('isclass',ec,'char'));
set(data.hQuiverc(i),'EdgeColor',col);

i = setdiff(1:length(data.hQuiverc),i);
set(data.hQuiverc(i),'FaceColor',col);

data.Options.Color = col;
set(data.hQuiverc,'UserData',data.Options);

mqSetColorControl(data);

guidata(obj,data);


% -------------------------------------------------
function mqSetShow(obj,eventdata)

data = guidata(obj);

shw = str2num(get(obj,'String'));
if (isempty(shw) | (shw < 0) | (shw > 100)),
    set(obj,'String',num2str(data.Options.Show*100));
else
    data.Options.Show = shw/100;
    data = mqUpdateVectors(data);
end;

guidata(obj,data);    

% -------------------------------------------------
function mqOffsetLength(obj,eventdata, type)

data = guidata(obj);
opts = data.Options;

switch type,
 case '+',
  opts.ScaleFactor = opts.ScaleFactor * 1.25;
 case '-',
  opts.ScaleFactor = opts.ScaleFactor * 0.8;
 case 'auto',
  opts.ScaleFactor = 1;
  opts.AbsScale = [];
end;

data.Options = opts;
data = mqUpdateVectors(data);

set(data.LengthEdit,'String',...
                  sprintf('%.3f',data.Options.ScaleFactor*data.Options.Scale));

guidata(obj,data);
  

% -------------------------------------------------
function mqSetLength(obj,eventdata)

data = guidata(obj);
opts = data.Options;

scale = str2num(get(data.LengthEdit,'String'));
if (~isempty(scale)),
    opts.AbsScale = scale;
    opts.ScaleFactor = 1;
end;

data.Options = opts;
data = mqUpdateVectors(data);
guidata(obj,data);



% -------------------------------------------------
function mqTruncate(obj,eventdata, type)

data = guidata(obj);

tt = get(data.TruncateTypeList,'Value');
if (tt == 1),
    trunc = data.Options.ScaleRange(2)*100;

    switch type,
     case '+',
      trunc = trunc + 5;
      if (trunc > 100),
          trunc = 100;
      end;
     case '-',
      trunc = trunc - 5;
      if (trunc < 0),
          trunc = 0;
      end;
     case 'set',
      t = str2num(get(obj,'String'));
      if (~isempty(t) & (t >= 0) & (t <= 100)),
          trunc = t;
      end;
     case 'off',
      trunc = 100;
    end;

    set(data.TruncateLengthEdit,'String',num2str(trunc));
    data.Options.ScaleRange(2) = trunc/100;
    data.Options.AbsScaleRange = [];
else
    trunc = data.Options.AbsScaleRange(2);
    maxlen = sqrt(max(data.Options.u(:).^2 + data.Options.v(:).^2));

    switch type,
     case '+',
      trunc = trunc * 1.25;
      if (trunc > maxlen),
          trunc = Inf;
      end;
     case '-',
      if (~isfinite(trunc)),
          trunc = maxlen*0.8;
      else
          trunc = trunc * 0.8;
      end;
     case 'set',
      t = str2num(get(obj,'String'));
      if (~isempty(t) & (t >= 0)),
          trunc = t;
      end;
     case 'off',
      trunc = Inf;
    end;

    set(data.TruncateLengthEdit,'String',num2str(trunc));
    data.Options.AbsScaleRange(2) = trunc;
    data.Options.ScaleRange = [];
end;

data = mqUpdateVectors(data);
set(data.LengthEdit,'String',...
                  num2str(data.Options.ScaleFactor*data.Options.Scale));

guidata(obj,data);

% -------------------------------------------------
function mqSetTruncateType(obj,eventdata)

data = guidata(obj);

len = sqrt(data.Options.u.^2 + data.Options.v.^2);

switch get(obj,'Value'),
 case 1,                                % relative scaling
  prc = sum(len(:) <= data.Options.AbsScaleRange(2));

  data.Options.ScaleRange = [0 diground(prc/prod(size(len)),0.01)];
  str = num2str(data.Options.ScaleRange(2)*100);

  data.Options.AbsScaleRange = [];


 case 2,                                % absolute scaling
  if (data.Options.ScaleRange(2) == 1),
      trunc = Inf;
  else
      trunc = prctile(len(:),data.Options.ScaleRange(2)*100);
  end;

  data.Options.AbsScaleRange = [0 trunc];
  str = num2str(trunc);

  data.Options.ScaleRange = [];
end;

set(data.TruncateLengthEdit,'String',str);
data = mqUpdateVectors(data);
guidata(obj,data);



% -------------------------------------------------
function mqToggleHeads(obj,eventdata)

data = guidata(obj);
data.Options.NoHeads = ~data.Options.NoHeads;

set(obj,'Value',~data.Options.NoHeads);

data = mqUpdateVectors(data);
guidata(obj,data);


% -------------------------------------------------
function mqSetHeadSize(obj,eventdata)

data = guidata(obj);

hs = get(obj,'Value');
data.Options.HeadSize = hs;

data = mqUpdateVectors(data);
guidata(obj,data);


% -------------------------------------------------
function mqSetHeadRange(obj,eventdata, elem,type)

data = guidata(obj);

edit = [data.HeadLoSizeEdit data.HeadHiSizeEdit];
tt = get(data.HeadRangeTypeList,'Value');

if (tt == 1),
    rng = data.Options.HeadRange*100;
    rngmax = [0 100];
    if (isempty(rng)),
        rng = rngmax;
    end;
else
    rng = data.Options.AbsHeadRange;
    rngmax = [0 Inf];
    if (isempty(rng)),
        rng = rngmax;
    end;
end;

switch type,
 case '+',
  if (tt == 1),
      rng(elem) = rng(elem) + 5;
  else
      rng(elem) = rng(elem) * 1.25;
  end;


 case '-',
  if (tt == 1),
      rng(elem) = rng(elem) - 5;
  else
      rng(elem) = rng(elem) * 0.8;
  end;

 case 'set',
  val = str2num(get(edit(elem),'String'));
  if (~isempty(val)),
      rng(elem) = val;
  end;
end;

if (rng(elem) < rngmax(1)),
    rng(elem) = rngmax(1);
elseif (rng(elem) > rngmax(2)),
    rng(elem) = rngmax(2);
end;

if (rng(1) > rng(2)),
    rng(3-elem) = rng(elem);
end;

set(edit(elem),'String',rng(elem));

if (tt == 1),
    data.Options.HeadRange = rng/100;
    data.Options.AbsHeadRange = [];
else
    data.Options.HeadRange = [];
    data.Options.AbsHeadRange = rng;
end;

data = mqUpdateVectors(data);
guidata(obj,data);

% -------------------------------------------------
function mqSetHeadRangeType(obj,eventdata)

data = guidata(obj);
len = flatten(sqrt(data.Options.u.^2 + data.Options.v.^2));
n = length(len);

switch get(obj,'Value'),
 case 1,                                % relative scaling
  lo = sum(len <= data.Options.AbsHeadRange(1))/n;
  hi = sum(len <= data.Options.AbsHeadRange(2))/n;

  data.Options.HeadRange = diground([lo hi],0.01);
  strlo = num2str(data.Options.HeadRange(1)*100);
  strhi = num2str(data.Options.HeadRange(2)*100);

  data.Options.AbsHeadRange = [];


 case 2,                                % absolute scaling
  rng = prctile(len,data.Options.HeadRange*100);
  if (data.Options.HeadRange(2) == 1),
      rng(2) = Inf;
  end;

  data.Options.AbsHeadRange = rng;
  strlo = num2str(rng(1));
  strhi = num2str(rng(2));

  data.Options.HeadRange = [];
end;

set(data.HeadLoSizeEdit,'String',strlo);
set(data.HeadHiSizeEdit,'String',strhi);
data = mqUpdateVectors(data);
guidata(obj,data);
    

% -------------------------------------------------
function mqResetHeadRange(obj,eventdata)

data = guidata(obj);
if (get(data.HeadRangeTypeList,'Value') == 1),
    data.Options.HeadRange = [0 1];
    data.Options.AbsHeadRange = [];
    set(data.HeadLoSizeEdit,'String','0');
    set(data.HeadHiSizeEdit,'String','1');
else
    data.Options.HeadRange = [];
    data.Options.AbsHeadRange = [0 Inf];
    set(data.HeadLoSizeEdit,'String','0');
    set(data.HeadHiSizeEdit,'String','Inf');
end;

data = mqUpdateVectors(data);
guidata(obj,data);

    



% -------------------------------------------------
function mqToggleTails(obj,eventdata)

data = guidata(obj);
data.Options.NoTails = ~data.Options.NoTails;

set(obj,'Value',~data.Options.NoTails);

data = mqUpdateVectors(data);
guidata(obj,data);


% -------------------------------------------------
function mqSetUICallbacks(data)

set(data.ColorList,'Callback',{@mqSetColor,'name'});
set(data.RGBEdit,'Callback',{@mqSetColor,'RGB'});
set(data.ShowEdit,'Callback',@mqSetShow);
set(data.LongerButton,'Callback',{@mqOffsetLength,'+'});
set(data.ShorterButton,'Callback',{@mqOffsetLength,'-'});
set(data.LengthEdit,'Callback',@mqSetLength);
set(data.ScaleAutoButton,'Callback',{@mqOffsetLength,'auto'});
set(data.TruncateUpButton,'Callback',{@mqTruncate,'+'});
set(data.TruncateDownButton,'Callback',{@mqTruncate,'-'});
set(data.TruncateLengthEdit,'Callback',{@mqTruncate,'set'});
set(data.TruncateTypeList,'Callback',@mqSetTruncateType);
set(data.TruncateOffButton,'Callback',{@mqTruncate,'off'});
set(data.HeadsOnCheck,'Callback',@mqToggleHeads);
set(data.HeadSizeSlider,'Callback',@mqSetHeadSize);
set(data.HeadLoUpButton,'Callback',{@mqSetHeadRange,1,'+'});
set(data.HeadLoDownButton,'Callback',{@mqSetHeadRange,1,'-'});
set(data.HeadLoSizeEdit,'Callback',{@mqSetHeadRange,1,'set'});
set(data.HeadHiUpButton,'Callback',{@mqSetHeadRange,2,'+'});
set(data.HeadHiDownButton,'Callback',{@mqSetHeadRange,2,'-'});
set(data.HeadHiSizeEdit,'Callback',{@mqSetHeadRange,2,'set'});
set(data.HeadRangeTypeList,'Callback',@mqSetHeadRangeType);
set(data.HeadRangeResetButton,'Callback',@mqResetHeadRange);
set(data.TailsOnCheck,'Callback',@mqToggleTails);


% -------------------------------------------------
function mqSetUIStart(data)

opts = data.Options;

mqSetColorControl(data);

set(data.LengthEdit,'String',num2str(opts.ScaleFactor*opts.Scale));
set(data.ShowEdit,'String',num2str(opts.Show*100));
if (~isempty(opts.AbsScaleRange)),
    set(data.TruncateLengthEdit,'String',num2str(opts.AbsScaleRange(2)));
    set(data.TruncateTypeList,'Value',2);
else,
    set(data.TruncateLengthEdit,'String',...
                      sprintf('%.1f',opts.ScaleRange(2)*100));
    set(data.TruncateTypeList,'Value',1);
end;

set(data.HeadsOnCheck,'Value',~opts.NoHeads);
set(data.HeadSizeSlider,'Value',opts.HeadSize);
if (~isempty(opts.AbsHeadRange)),
    set(data.HeadLoSizeEdit,'String',num2str(opts.AbsHeadRange(1)));
    set(data.HeadHiSizeEdit,'String',num2str(opts.AbsHeadRange(2)));
    set(data.HeadRangeTypeList,'Value',2);
else
    set(data.HeadLoSizeEdit,'String',num2str(opts.HeadRange(1)*100));
    set(data.HeadHiSizeEdit,'String',num2str(opts.HeadRange(2)*100));
    set(data.HeadRangeTypeList,'Value',1);
end;

set(data.TailsOnCheck,'Value',~opts.NoTails);

% -------------------------------------------------
function mqSetColorControl(data)

RGB = [0 0 1; 0 1 0; 1 0 0; 0 1 1; 1 0 1; 1 1 0; 1 1 1; 0 0 0];
Names = {'Blue'; 'Green'; 'Red'; 'Cyan'; 'Magenta'; 'Yellow'; ...
         'White'; 'Black'};

switch get(data.hQuiverc(1),'Type'),
 case 'line',
  col = get(data.hQuiverc(1),'Color');
 case 'patch',
  ecol = get(data.hQuiverc,'EdgeColor');
  k = find(~cellfun('isclass',ecol,'char'));
  if (~isempty(k)),
      col = ecol{k(1)};
  else
      fcol = get(data.hQuiverc,'FaceColor');
      k = find(~cellfun('isclass',ecol,'char'));
      if (isempty(k)),
          return;
      end;
      col = fcol{k(1)};
  end;
end;

ind = find(all(RGB == repmat(col,[size(RGB,1) 1]),2));

if (~isempty(ind)),
    set(data.ColorList,'Value',ind);
else
    set(data.ColorList,'Value',length(Names)+1);
end;
str = sprintf('[%.1f %.1f %.1f]',col);

set(data.RGBEdit,'String',str);


    