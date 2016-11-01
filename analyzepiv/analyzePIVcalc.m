function [data,calcguifcns,calcfcns] = analyzePIVcalc(data)

calcguifcns.apAddCalc = @apAddCalc;
calcguifcns.apScrollCalc = @apScrollCalc;
calcguifcns.apSlideCalcScroll = @apSlideCalcScroll;
calcguifcns.apSetRgnNames = @apSetRgnNames;
calcguifcns.apUpdateCalc = @apUpdateCalc;
calcguifcns.apScrollCalc = @apScrollCalc;
calcguifcns.apDeleteCalc = @apDeleteCalc;
calcguifcns.apCalculate = @apCalculate;
calcguifcns.apSaveData = @apSaveData;

% second column defines what parameters the functions need
% 1 = boundary only
% 2 = internal vectors only
% 3 = both
calcfcns = {@apCalcCirc 1 'circ'; ...
            @apCalcStepCirc 3 'scirc'; ...
            @apCalcMeanMag 2 'magmean'; ...
            @apCalcMaxMag 2 'magmax'; ...
            @apCalcMeanUVW 2 'uvmean'; ...
            @apCalcMeanAng 2 'angmean'; ...
            @apCalcFlux 1 'flux'};

set(data.RgnDataScroll,'Callback',@apSlideCalcScroll);

% setup the calculation box rows
gap = 3;

pos = get(data.RgnName0,'Position');
ext = get(data.RgnFrame,'Position');

h = pos(4);
% we use the center because some of the different controls
% are different heights - this allows us to keep them centered
% w/r/t each other
startpos = pos(2) - h/2 - gap;
newrowctr = startpos;

row = 1;
while (newrowctr - h/2 > ext(2) + gap),
    % add another row of controls for this region
    GUI.Num = copyobj(data.RgnNum0,data.Panel);
    GUI.Name = copyobj(data.RgnName0,data.Panel);
    GUI.Fcns(1) = copyobj(data.DataCalcCirc0,data.Panel);
    GUI.Fcns(2) = copyobj(data.DataCalcSCirc0,data.Panel);
    GUI.Fcns(3) = copyobj(data.DataCalcMag0,data.Panel);
    GUI.Fcns(4) = copyobj(data.DataCalcMaxMag0,data.Panel);
    GUI.Fcns(5) = copyobj(data.DataCalcUV0,data.Panel);
    GUI.Fcns(6) = copyobj(data.DataCalcAng0,data.Panel);
    GUI.Fcns(7) = copyobj(data.DataCalcMomFlux0,data.Panel);
    
    GUI.Properties = copyobj(data.RgnPropertiesButton0,data.Panel);

    % put the new row in the right position
    pos = get(GUI.Num,'Position');
    pos(2) = newrowctr - pos(4)/2;
    set(GUI.Num, 'Position',pos);

    pos = get(GUI.Name,'Position');
    pos(2) = newrowctr - pos(4)/2;
    set(GUI.Name,'Callback',@apRgnNameEdit,...
                 'Position',pos);

    for i = 1:length(calcfcns),
        pos = get(GUI.Fcns(i),'Position');
        pos(2) = newrowctr - pos(4)/2;
        set(GUI.Fcns(i),'Callback',@apDataCalcFcnClick,...
                        'Position',pos);
    end;
    pos = get(GUI.Properties,'Position');
    pos(2) = newrowctr - pos(4)/2;
    set(GUI.Properties, 'Callback',@apDataRgnPropClick,...
                      'Position',pos, 'Visible','off');

    data.Calc.GUI(row) = GUI;
    row = row+1;

    newrowctr = startpos - (row-1)*(h + gap);
end;
data.Calc.maxVisibleRows = row-1;
data.Calc.TopRow = 1;
data.Calc.doFcns = zeros(0,length(calcfcns));

% -------------------------------------------------
function data = apAddCalc(data)
% now somewhat synonymous with adding a region - regions and
% calculations used to be separated, but that turned out to make
% the GUI really complicated

% add the new elements to specify the calculations for the new
% regions
newfcns = zeros(size(data.Regions,1)-size(data.Calc.doFcns,1), ...
                length(data.calcFcns));
newfcns(1) = 1;     % circulation calculation is on by default
data.Calc.doFcns = cat(1,data.Calc.doFcns,newfcns);

data = apScrollCalc(data,'visible',size(data.Regions,1));


% -------------------------------------------------
function data = apUpdateCalc(data)

for i = 1:data.Calc.maxVisibleRows,
    rgn = data.Calc.TopRow + i - 1;

    if (rgn <= size(data.Regions,1)),
        set(data.Calc.GUI(i).Num,'String',strcat(num2str(rgn),'.'),...
                          'Visible','on','UserData',rgn);
        set(data.Calc.GUI(i).Name,'String',data.rgnNames{rgn}, ...
                          'Visible','on','UserData',rgn);
        set(data.Calc.GUI(i).Properties, 'Visible','on','UserData',rgn);
        
        % some region types can't do all of the functions
        switch (data.Regions(rgn,data.curFrame).type),
         case {'line','poly'},
          enab = {'on', 'off', 'on', 'on', 'on', 'on','on'};
         otherwise,
          enab = {'on', 'on', 'on', 'on', 'on', 'on','on'};
        end;

        for j = 1:length(data.calcFcns),
            set(data.Calc.GUI(i).Fcns(j),...
                'Value', data.Calc.doFcns(rgn,j),...
                'Visible','on','Enable',enab{j},...
                'UserData',[rgn j]);
        end;
    else
        h = [data.Calc.GUI(i).Num data.Calc.GUI(i).Name ...
             data.Calc.GUI(i).Fcns data.Calc.GUI(i).Properties];
        set(h, 'Visible','off');
    end;
end;

% check if we need to enable the slider
if (data.Calc.maxVisibleRows < size(data.Regions,1)),
    % slider position determines the *top* row in the calculation
    % list.  So the max position has to be less than the maximum
    % number of regions.
    maxscroll = size(data.Regions,1) - data.Calc.maxVisibleRows + 1;
    ss = 1/(maxscroll-1);
    ss(2) = max(ss,0.1);           % 10% or 1 element, whichever is bigger
    val = get(data.RgnDataScroll,'Value');
    if (val > maxscroll),
        val = maxscroll;
    end;
    set(data.RgnDataScroll,'Visible','on', ...
                      'Value',val,...
                      'Min',1, 'Max', maxscroll, ...
                      'SliderStep',ss);
else
    set(data.RgnDataScroll,'Visible','off');
end;


% -------------------------------------------------
function data = apScrollCalc(data,func,d)

switch lower(func),
 % scroll until a particular item is visible
 case 'visible',
  % check if it's off screen and set the top row appropriately
  if (d < data.Calc.TopRow),
      top = d;
  elseif (d > data.Calc.TopRow + data.Calc.maxVisibleRows - 1),
      top = d - data.Calc.maxVisibleRows + 1;
  else
      top = data.Calc.TopRow;
  end;

 % set the top item
 case 'top',
  top = d;
end;

% check if the top index is too high
if (top + data.Calc.maxVisibleRows - 1 > size(data.Regions,1)),
    % set top to be the first row that means that all visible rows
    % will be filled, or 1, whichever is higher (in order to deal
    % with the times we call apScrollCalc when there actually aren't
    % enough elements to scroll)
    top = max(1,size(data.Regions,1) - data.Calc.maxVisibleRows + 1);
end;

data.Calc.TopRow = top;
data = apUpdateCalc(data);

% set the slider to be in this position
% or at position 1, if we don't have enough elements to scroll
maxscroll = get(data.RgnDataScroll,'Max');
set(data.RgnDataScroll,'Value',max(maxscroll - top + 1, 1));

% -------------------------------------------------
function apSlideCalcScroll(obj,eventdata)

data = guidata(obj);

% NB: Vertical sliders scroll (for some stupid reason) from the bottom
% so we have to invert the value
top = get(obj,'Value');
maxscroll = get(data.RgnDataScroll,'Max');

data = apScrollCalc(data, 'top',round(maxscroll - top + 1));

guidata(obj,data);

% -------------------------------------------------
function apDataCalcFcnClick(obj,eventdata)

data = guidata(obj);

dofcn = get(obj,'Value');
ind = get(obj,'UserData');

data.Calc.doFcns(ind(1),ind(2)) = dofcn;

guidata(obj,data);

% -------------------------------------------------
function apRgnNameEdit(obj,eventdata)

data = guidata(obj);

rgn = get(obj,'UserData');

data.rgnNames{rgn} = get(obj,'String');
set(data.Regions(rgn,data.curFrame).label, 'String', data.rgnNames{rgn});

guidata(obj,data);

% -------------------------------------------------
function data = apDeleteCalc(data,rgn)
% really means the region was deleted

data.Calc.doFcns = data.Calc.doFcns([1:rgn-1 rgn+1:end],:);

% we use scroll here instead of update because we want to scroll
% the box up if we're at the end of the region list, but we still
% have more regions than rows in the box.  This prevents us from
% displaying a blank row at the bottom
data = apScrollCalc(data,'top',data.Calc.TopRow);

% -------------------------------------------------
function apDataRgnPropClick(obj,eventdata)

data = guidata(obj);

rgnnum = get(obj,'UserData');

% find when the region is defined
isdef = find(cat(2,data.Regions(rgnnum,:).status) > 0);

if (data.Regions(rgnnum,data.curFrame).status == 0),
    % find the closest frame in which the region is defined
    [m,fr] = min(abs(isdef-data.curFrame));
    fr = isdef(fr);

    guidata(obj,data);
    feval(data.generalFcns.apSetFrame, obj,eventdata, 'set',fr);
    data = guidata(obj);
end;

if (data.Regions(rgnnum,data.curFrame).status == 2),
    data = feval(data.rgnFcns.apSetCurRgn, data, rgnnum);
end;

dlg = openfig('analyzePIVproperties.fig', 'new');
dlgdata = guihandles(dlg);

set(dlgdata.RgnNameEdit, 'String',data.rgnNames{rgnnum});
set(dlgdata.StartFrameEdit, 'String',num2str(isdef(1)));
set(dlgdata.EndFrameEdit, 'String',num2str(isdef(end)));
set(dlgdata.ToBeginningButton, 'Callback', @apRgnPropsToBeginning);
set(dlgdata.ToEndButton, 'Callback', @apRgnPropsToEnd);
set(dlgdata.CancelButton,'Callback','delete(get(gcbo,''Parent''));');
set(dlgdata.OKButton,'Callback','uiresume;');
dlgdata.nFrames = data.nFrames;
guidata(dlg,dlgdata);

uiwait(dlg);

if (ishandle(dlg)),
    rgnname = get(dlgdata.RgnNameEdit, 'String');
    if (~strcmp(rgnname,data.rgnNames{rgnnum})),
        data.rgnNames{rgnnum} = rgnname;
        data = apUpdateCalc(data);
    end;

    % reset the start frame
    rgnstart = str2num(get(dlgdata.StartFrameEdit, 'String'));
    rgnstart = max(rgnstart,1);
    if (rgnstart < isdef(1)),
        for i = rgnstart:isdef(1)-1,
            data.Regions(rgnnum,i).status = 1;
        end;
    end;

    % reset the end frame
    rgnend = str2num(get(dlgdata.EndFrameEdit, 'String'));
    rgnend = min(rgnend,data.nFrames);
    if (rgnend > isdef(end)),
        for i = isdef(end)+1:rgnend,
            data.Regions(rgnnum,i).status = 1;
        end;
    end;
    
    data = feval(data.rgnFcns.apInterpRgn,data);
    delete(dlg);
end;
    
guidata(obj,data);    


% -------------------------------------------------
function apRgnPropsToBeginning(obj,eventdata)

dlgdata = guidata(obj);

set(dlgdata.StartFrameEdit,'String','1');
guidata(obj,dlgdata);

% -------------------------------------------------
function apRgnPropsToEnd(obj,eventdata)

dlgdata = guidata(obj);

set(dlgdata.EndFrameEdit,'String',num2str(dlgdata.nFrames));
guidata(obj,dlgdata);



% -------------------------------------------------
function [data, outsimple, out, err, units, N, ctr, name] = ...
    apCalculate(obj, eventdata, calcind)
% actually do the calculations
% originally set up to do individual calculations, but now
% almost entirely used for recalculating for all regions at once

if (~ishandle(obj)),
    data = obj;
else
    data = guidata(obj);
end;

bStereo = data.PIV.isStereo;

if (isempty(calcind)),
    if (ishandle(obj)),
        calcind = get(obj,'UserData');
    else
        calcind = 1:size(data.Regions,1);
    end;
end;

if (size(data.PIV.x,3) > 1),
    px = data.PIV.x(:,:,data.curFrame);
    py = data.PIV.y(:,:,data.curFrame);
else
    px = data.PIV.x(:,:,1);
    py = data.PIV.y(:,:,1);
end;
if (isfield(data.PIV,'us'))
    pu = data.PIV.us(:,:,data.curFrame) - data.subtractVector(1);
    pv = data.PIV.vs(:,:,data.curFrame) - data.subtractVector(2);
else
    pu = data.PIV.u(:,:,data.curFrame) - data.subtractVector(1);
    pv = data.PIV.v(:,:,data.curFrame) - data.subtractVector(2);
end

if bStereo
    if isfield(data.PIV,'ws')
        pw = data.PIV.ws(:,:,data.curFrame) - data.subtractVector(3);
    else
        pw = data.PIV.w(:,:,data.curFrame) - data.subtractVector(3);
    end
else
    pw = [];
end

d = nanmean(abs([px(1,2)-px(1) py(2,1)-py(1)]));

% set things empty at the beginning
out = {};
outsimple = {};
err = {};
units = {};
N = {};
ctr = {};
name = {};

% run through all the regions
calcnum = 1;                        % tracks which calculations are checked
for i = 1:length(calcind),
    rgn = calcind(i);

    % and run through the possible things we can calculate for
    % each region
    for fcn = 1:length(data.calcFcns),
        docalc = data.Calc.doFcns(rgn,fcn);
        frm = data.curFrame;

        % get the contour
        x0 = data.Regions(rgn,frm).x;
        y0 = data.Regions(rgn,frm).y;
        
        if (docalc & (data.Regions(rgn,frm).status > 0)),
            % check what type of data we need
            % 1 or 3 means velocities along the contour
            % also, a line has no interior data, so interpolate its
            % contour regardless
            if ((data.calcFcns{fcn,2} == 1) | (data.calcFcns{fcn,2} == 3) | ...
                strcmp(data.Regions(rgn,frm).type,'line')),
                s0 = [0 cumsum(sqrt(diff(x0).^2 + diff(y0).^2))];

                s1 = 0:d:s0(end);
                xd = interp1(s0,x0, s1);
                yd = interp1(s0,y0, s1);

                ud = interp2(px,py,pu, xd,yd);
                vd = interp2(px,py,pv, xd,yd);
                if bStereo, 
                    wd = interp2(px,py,pw, xd, yd); 
                else, 
                    wd = []; 
                end
            else
                xd = x0;
                yd = y0;
                ud = [];
                vd = [];
                wd = [];
            end;

            % for a line, set the interior data to be the data along
            % the contour
            if (strcmp(data.Regions(rgn,frm).type,'line')),
                xv = xd;
                yv = yd;
                uv = ud;
                vv = vd;
                if (bStereo),
                    wv = wd;
                else
                    wv = [];
                end;
            elseif ((data.calcFcns{fcn,2} == 2) | ...
                    (data.calcFcns{fcn,2} == 3))
                % otherwise, grab the interior vectors
                [xv,yv,uv,vv,wv] = apGetVecInside(data,x0,y0);
            else
                xv = [];
                yv = [];
                uv = [];
                vv = [];
                wv = [];
            end;

            % and then actually run the calculation
            [out1,err1,units1,N1,ctr1] = ...
                feval(data.calcFcns{fcn,1},data,rgn,xd,yd,ud,vd,wd, ...
                      xv,yv,uv,vv,wv, data.Units);
        else
            out1 = [];
            err1 = [];
            units1 = '';
            out2 = [];
            N1 = 0;
            ctr1 = [];
        end;

        if (docalc),
            % some calculation functions return more sophisticated output
            % for csv output, we need to get the most important bit, which
            % is the first element in whatever it is
            if (iscell(out1)),
                outsimple{calcnum,1} = out1{1};
            elseif (isstruct(out1)),
                out2 = struct2cell(out1);
                outsimple{calcnum,1} = out2{1};
            else
                outsimple{calcnum,1} = out1;
            end;
            
            out{calcnum,1} = out1;
            err{calcnum,1} = err1;
            units{calcnum,1} = units1;
            N{calcnum,1} = N1;
            ctr{calcnum,1} = ctr1;
            
            % build up a name for this calculation by stringing together
            % the region name and the function name
            name{calcnum,1} = strcat(data.rgnNames{i}, ...
                                     data.calcFcns{fcn,3});

            % remove any spaces from the name
            k = find(name{calcnum,1} ~= ' ');
            name{calcnum,1} = name{calcnum,1}(k);

            calcnum = calcnum+1;
        end;
    end;
end;

if (ishandle(obj)),
    guidata(obj,data);
end;



% -------------------------------------------------
function [doclose,data] = apSaveData(data)

curFr = data.curFrame;

if (all(data.Calc.doFcns == 0)),
    doclose = 1;
    return;
end;

if (isempty(data.SaveDataFile)),
    filter = {'*.csv','Comma separated values for Excel (*.csv)';
              '*.mat','Matlab data file (*.mat)'};
    filtval = {'csv','matlab'};

    [file,path,filtind] = uiputfile(filter, 'Save data to file...');

    if (filtind ~= 0),
        file = fullfile(path,file);
        data.SaveDataFile = file;
        data.SaveDataType = filtval{filtind};
    else
        doclose = 1;
        data.SaveDataFile = '';
        data.SaveDataType = [];
        return;
    end;
end;

% interpolate regions between keyframes before we calculate
data = feval(data.rgnFcns.apInterpRgn,data);

try,
    timedWaitBar(0,'Calculating...');
    for i = 1:data.nFrames,
        data.curFrame = i;
        [data,val00,valfull00,err00,units0,N00,ctr00,names0] = ...
            apCalculate(data,[],[]);

        if (strcmp(data.SaveDataType,'csv')),
            val0(:,i) = val00;
        else
            val0(:,i) = valfull00;
        end;
        err0(:,i) = err00;
        N0(:,i) = N00;
        ctr0(:,i) = ctr00;

        if (~timedWaitBar(i/data.nFrames)),
            doclose = 0;
            return;
        end;
    end;
    timedWaitBar(1);
catch,
    fprintf('Error during calculation.  Attempting to recover...');
    timedWaitBar(1);
end;

switch (data.SaveDataType),

%------------------------------------------------
% Full Matlab file output - more data, but only readable in
% Matlab
 case 'matlab',
  if (exist(data.SaveDataFile)),
      existing = who('-file',data.SaveDataFile);
      opt = {'-append'};
  else
      existing = {};
      opt = {};
  end;

  names = {};
  overwrite = {};
  val = {};
  err = {};
  errnames = {};
  for i = 1:length(names0),
      if (~isvarname(names0{i})),       % check if it's a valid var name
          oldname = names0{i};

          % put a V at the beginning of the name if we start with a
          % digit
          names0{i} = regexprep(names0{i},'^(\d+.*)$','V$1','tokenize');
          % replace unsavory characters with underscores
          names0{i} = regexprep(names0{i},'[-+~!@#$%^&*]','_');
          if (length(names0{i}) > namelengthmax),
              names0{i} = names0{i}(1:namelengthmax);
          end;

          warning('Variable name %s is invalid.  Changed to %s.',...
                  oldname,names0{i});
      end;

      names{end+1} = names0{i};
      units.(names0{i}) = units0{i};
      val{end+1} = val0(i,:);
      if (any(~cellfun('isempty',err0(i,:)))),
          names{end+1} = [names0{i} 'err'];
          val{end+1} = err0(i,:);
      end;
      names{end+1} = [names0{i} 'N'];
      val{end+1} = N0(i,:);

      names{end+1} = [names0{i} 'ctr'];
      val{end+1} = ctr0(i,:);
  end;

  names{end+1} = 'Units';
  val{end+1} = units;

  for i = 1:length(names),
      if (~isempty(strmatch(names{i},existing,'exact'))),
          overwrite{end+1} = names{i};
      end;
  end;

  if (~isempty(overwrite)),
      over = sprintf('%s ',overwrite{:});
      prompt = sprintf('Overwrite %d variables?\n(%s)',...
                       length(overwrite),over(1:end-1));
      butt = questdlg(prompt, 'Overwrite?', ...
                      'Yes','No and Continue','No and Quit',...
                      'No and Continue');
      switch butt,
       case 'No and Continue',
        doclose = 0;
        return;
       case 'No and Quit',
        doclose = 1;
        return;
      end;
  end;

  for iYind = 1:size(names,2)
      nmind{iYind} = iYind;
  end

  for iCnt = 1:length(names)
      fcn = sprintf('%s=val{%d};', names{iCnt}, nmind{iCnt});
      eval(fcn);
  end
  for iCnt = 1:length(errnames)
      fcn = sprintf('%s=err{%d};', errnames{iCnt}, nmind{iCnt});
      eval(fcn);
  end
  
  
  save(data.SaveDataFile,names{:},errnames{:},opt{:});


%-----------------------------------------------------------
% comma separated value file output - readable by Excel
% and standard stats packages, but doesn't include all
% the information we can save in a Matlab file
 case 'csv',
  fid = fopen(data.SaveDataFile,'w');
  if (isfield(data.PIV,'files')),
      if (iscellstr(data.PIV.files)),
          fn1 = data.PIV.files{1};
          fn2 = data.PIV.files{end};
          fprintf(fid,'%% Data from %s files, %s - %s\n',data.PIV.fileType,...
                  fn1,fn2);
      else
          fprintf(fid,'%% Data from %s file %s\n',data.PIV.fileType,...
                  data.PIV.files);
      end;
  else
      fprintf(fid,'%% Data from Matlab base workspace\n');
  end;
  fprintf(fid,'%% %d x %d vectors, %d frames\n',size(data.PIV.u,1),...
          size(data.PIV.u,2), size(data.PIV.u,3));
  fprintf(fid,'%% %d regions\n',size(data.Regions,1));
  fprintf(fid,'%% ');
  for i = 1:size(data.Regions,1)-1,
      fprintf(fid,'Rgn %d = %s, ',i,data.Regions(i,1).type);
  end;
  fprintf(fid,'Rgn %d = %s\n',size(data.Regions,1),...
          data.Regions(end,1).type);
  fprintf(fid,'%% %s\n',datestr(clock,0));
  fprintf(fid,'\n\n');

  waserror = 0;
  try,
      head = {'Frame'};
      out = data.PIV.frames';
      units = {'fr'};
      col = 2;
      a = 2;
      for i = 1:length(names0),
          nm = names0{i};

          % some calculations return multiple values (mostly the mean UV)
          % we'll need more than one column for those
          ncol1 = max(cellfun('prodofsize',val0(i,:)));

          % set up our columns
          v1 = repmat(NaN,[length(data.PIV.frames) ncol1]);

          % only assign where the value isn't empty - otherwise we'll
          % get an error
          good = find(~cellfun('isempty',val0(i,:)));
          v1(good,:) = cat(1,val0{i,good});
          out(:,end+(1:ncol1)) = v1;

          if (ncol1 == 1),
              head{1,end+1} = nm;
              units(1,end+1) = units0(i,1);
          else
              for j = 1:ncol1,
                  head{1,end+1} = strcat(nm,num2str(j));
                  units(1,end+1) = units0(i,1);
              end;
          end;

          % attach columns with the error estimate, if we have it
          % again, if we had multiple columns of data, we could have multiple
          % columns of error
          if (any(~cellfun('isempty',err0(i,:)))),
              err1 = repmat(NaN,[size(out,1) ncol1]);
              err1(good,:) = cat(1,err0{i,:});
              
              % add "err" to the end of the name
              out(:,end+(1:ncol1)) = err1;
              if (ncol1 == 1),
                  head{1,end+1} = strcat(nm,'err');
                  units(1,end+1) = units0(i,1);
              else
                  for j = 1:ncol1,
                      head{1,end+1} = strcat(nm,'err',num2str(j));
                      units(1,end+1) = units0(i,1);
                  end;
              end;
          end;
          
          % set up a column for the number of vectors used in the
          % calculation
          N1 = cat(1,N0{i,:});
          out(:,end+1) = N1;
          head{1,end+1} = strcat(nm,'N');
          units{1,end+1} = '';

          % and also two columns for the x and y center of the region
          ctr1 = repmat(NaN,[size(out,1) 2]);
          ctr1(good,:) = cat(1,ctr0{i,:});
          out(:,end+(1:2)) = ctr1;
          head(1,end+(1:2)) = {strcat(nm,'ctrX') strcat(nm,'ctrY')};
          units(1,end+(1:2)) = cellstr(repmat(data.Units.pos,[2 1]))';

          % check for my own stupidity - if for some reason we ended up
          % with more columns in the output matrix than in the header,
          % then give a warning and generate something to fill in the
          % space
          if (size(head,2) ~= size(out,2)),
              warning('Internal error in column names.');
              if (size(head,2) > size(out,2)),
                  head = head(:,1:size(out,2));
              else
                  head(:,size(head,2)+1:size(out,2)) = ...
                      cellstr(repmat('?',[size(out,2)-size(head,2) 1]))';
              end;
          end;
          % same thing with the units variable
          if (size(units,2) ~= size(out,2)),
              warning('Internal error in column names.');
              if (size(units,2) > size(out,2)),
                  units = units(:,1:size(out,2));
              else
                  units(:,size(units,2)+1:size(out,2)) = ...
                      cellstr(repmat('?',[size(out,2)-size(units,2) 1]))';
              end;
          end;
      end;

  catch,
      % more checks for my own stupidity - we shouldn't be getting errors
      % while we do the calculations, but if we do, try to salvage the
      % file so that all the data aren't lost
      fprintf(['Error while calculating data file.  Attempting to ' ...
               'recover...\n']);
      
      if (exist('head') & (length(head) ~= size(out,2))),
          warning('Column names are messed up...');
      elseif (~exist('head')),
          warning('No column names...');
          head = {};
      end;
      if (exist('units') & (length(units) ~= size(out,2))),
          warning('Units are messed up...');
      elseif (~exist('units')),
          warning('No units...');
          units = {};
      end;
      
      waserror = 1;
  end;
  
  col = size(out,2);

  tplt = strcat(repmat('%s, ',[1 col-1]),' %s\n');
  fprintf(fid,tplt,head{:});
  tplt = strcat(repmat('[%s], ',[1 col-1]),' [%s]\n');
  fprintf(fid,tplt,units{:});

  tplt = strcat(repmat('%g, ',[1 col-1]),' %g\n');
  fprintf(fid,tplt,out');

  fclose(fid);
  
  if (waserror),
      rethrow(lasterror);
  end;
end;

data.curFrame = curFr;
doclose = 1;


% -------------------------------------------------
function [circ,err,units,N,ctr] = apCalcCirc(data,rgn,xd,yd,ud,vd,wd, ...
                                   xv,yv,uv,vv,ww, vecunits)
% wd and ww are not used because circulation is calculated in the xy plane.

% get the scale factor between position units and velocity length unit
[posscalefac,unitserror] = feval(data.generalFcns.apConvertUnits, ...
                                 vecunits.pos,vecunits.vel{1});

if (~unitserror),
    vecunits.pos = vecunits.vel{1};
end;

xd = xd.*posscalefac;
yd = yd.*posscalefac;

s = [0 cumsum(sqrt(diff(xd).^2 + diff(yd).^2))];

if length(xd) >= 3
    dxds = deriv(s,xd);
    dyds = deriv(s,yd);
    mag = sqrt(dxds.^2 + dyds.^2);
    dxds = dxds./mag;
    dyds = dyds./mag;

    circ = trapz(s, ud.*dxds + vd.*dyds);
    N = length(s);
else
    circ = NaN;
    N = length(xd);
end

% convert ctr back to the original position units
ctr = [nanmean(xd(:)) nanmean(yd(:))]/posscalefac;

units = '';
if (~isempty(vecunits.pos) & ~isempty(vecunits.vel)),
    if (strcmp(vecunits.pos,vecunits.vel{1})),
        units = sprintf('%s^2/%s',vecunits.vel{1},vecunits.vel{2});
    else
        units = sprintf('%s.%s/%s',vecunits.pos,vecunits.vel{1},...
                        vecunits.vel{2});
    end;
end;
err = [];

% -------------------------------------------------
function [scirc,err,units,N,ctr] = apCalcStepCirc(data,rgn,xd,yd,ud,vd,wd, ...
                                   xv,yv,uv,vv,wv, vecunits)
% wd and ww are not used because circulation is calculated in the xy plane.
if (size(data.PIV.x,3) > 1),
    px = data.PIV.x(:,:,data.curFrame);
    py = data.PIV.y(:,:,data.curFrame);
else
    px = data.PIV.x(:,:,1);
    py = data.PIV.y(:,:,1);
end;
if (isfield(data.PIV,'us'))    
    pu = data.PIV.us(:,:,data.curFrame);
    pv = data.PIV.vs(:,:,data.curFrame);
else
    pu = data.PIV.u(:,:,data.curFrame);
    pv = data.PIV.v(:,:,data.curFrame);
end

d = nanmean([flatten(diff(px,[],2)); flatten(diff(py,[],2))]);

x = data.Regions(rgn,data.curFrame).x;
y = data.Regions(rgn,data.curFrame).y;

k = convhull(x,y);                      % for dealing with weird polygons
ctrx = mean(x(k));
ctry = mean(y(k));

x = x-ctrx;
y = y-ctry;

rmax = max(sqrt(x.^2 + y.^2));
r = 0:d:rmax;
r = (r + (rmax-r(end)));
step = r/rmax;

[posscalefac,unitserror] = feval(data.generalFcns.apConvertUnits, ...
                                 vecunits.pos,vecunits.vel{1});

if (~unitserror),
    vecunits.pos = vecunits.vel{1};
end;

xd = {};
yd = {};
for i = 1:length(step),
    x0 = x*step(i) + ctrx;
    y0 = y*step(i) + ctry;

    s0 = [0 cumsum(sqrt(diff(x0).^2 + diff(y0).^2))];

    s1 = 0:d:s0(end);
    if (length(s1) >= 4),
        xd1 = interp1(s0,x0, s1);
        yd1 = interp1(s0,y0, s1);

        ud = interp2(px,py,pu, xd1,yd1);
        vd = interp2(px,py,pv, xd1,yd1);

        dxds = deriv(s1,xd1);
        dyds = deriv(s1,yd1);
        mag = sqrt(dxds.^2 + dyds.^2);
        dxds = dxds./mag;
        dyds = dyds./mag;

        circ(i) = trapz(posscalefac*s1, ud.*dxds + vd.*dyds);
        xd{i} = xd1;
        yd{i} = yd1;
    else
        circ(i) = NaN;
        xd{i} = [];
        yd{i} = [];
    end;
end;

[cmax,ind] = max(abs(circ));
scirc.circ = circ(ind);
scirc.circSteps = circ;
scirc.r = r;
scirc.boundx = xd;
scirc.boundy = yd;

N = length(xd{ind});
ctr = [ctrx ctry];

units = '';
if (~isempty(vecunits.pos) & ~isempty(vecunits.vel)),
    if (strcmp(vecunits.pos,vecunits.vel{1})),
        units = sprintf('%s^2/%s',vecunits.vel{1},vecunits.vel{2});
    else
        units = sprintf('%s.%s/%s',vecunits.pos,vecunits.vel{1},...
                        vecunits.vel{2});
    end;
end;
err = [];


% -------------------------------------------------
function [mn,err,units,N,ctr] = apCalcMeanUVW(data,rgn,xd,yd,ud,vd,wd, ...
                                   xv,yv,uv,vv,wv, vecunits)
% should work with stereo and non-stereo results.
if data.PIV.isStereo
    N = sum(isfinite(uv) & isfinite(vv) & isfinite(wv));
else
    N = sum(isfinite(uv) & isfinite(vv));
end
mn = [nanmean(uv(:)) nanmean(vv(:)) nanmean(wv(:))];
if (N > 1),
    err = [nanstd(uv(:)) nanstd(vv(:)) nanstd(wv(:))] ./ ...
          sqrt(sum(isfinite(uv(:))));
else
    err = [NaN NaN NaN];
end;
if (~isempty(vecunits.vel)),
    units = [vecunits.vel{1} '/' vecunits.vel{2}];
else
    units = '';
end;
ctr = [nanmean(xd(:)) nanmean(yd(:))];



% -------------------------------------------------
function [mn,err,units,N,ctr] = apCalcMeanMag(data,rgn,xd,yd,ud,vd,wd, ...
                                   xv,yv,uv,vv,wv, vecunits)
if ~data.PIV.isStereo
    mag = sqrt(uv(:).^2 + vv(:).^2);
else
    mag = sqrt(uv(:).^2 + vv(:).^2 + wv(:).^2);
end
mn = nanmean(mag);
err = [nanstd(mag)]./sqrt(sum(isfinite(mag)));

if (~isempty(vecunits.vel)),
    units = [vecunits.vel{1} '/' vecunits.vel{2}];
else
    units = '';
end;

if ~data.PIV.isStereo
    N = sum(isfinite(uv) & isfinite(vv));
else
    N = sum(isfinite(uv) & isfinite(vv) & isfinite(wv));
end

ctr = [nanmean(xd(:)) nanmean(yd(:))];

% -------------------------------------------------
function [mn,err,units,N,ctr] = apCalcMeanAng(data,rgn,xd,yd,ud,vd,wd, ...
                                        xv,yv,uv,vv,wv, vecunits)
% doesn't work on the Z component.
% Mean of angle data from Fisher (1993)
ang = atan2(vv(:),uv(:));
C = nanmean(cos(ang));
S = nanmean(sin(ang));
r = sqrt(C.^2 + S.^2);
mn = atan2(S,C)*180/pi;

% Circular standard error Fisher (1993) eq 4.42
N = sum(isfinite(uv) & isfinite(vv));

if (N > 2),
    kappa = angkappa(ang);
    err = 1/sqrt(sum(isfinite(ang(:)))*r*kappa) * 180/pi;
else
    err = NaN;
end;
units = 'deg';
ctr = [nanmean(xd(:)) nanmean(yd(:))];

% -------------------------------------------------
function [mn,err,units,N,ctr] = apCalcMed(data,rgn,xd,yd,ud,vd,wd, ...
                                   xv,yv,uv,vv,wv, vecunits)

mn = [nanmedian(uv(:)) nanmedian(vv(:)) nanmedian(wv(:))];
err = [];
if (~isempty(vecunits.vel)),
    units = [vecunits.vel{1} '/' vecunits.vel{2}];
else
    units = '';
end;
N = sum(isfinite(uv) & isfinite(vv));
ctr = [nanmean(xd(:)) nanmean(yd(:))];

% -------------------------------------------------
function [mn,err,units,N,ctr] = apCalcMin(data,rgn,xd,yd,ud,vd,wd, ...
                                   xv,yv,uv,vv,wv, vecunits)

mn = [nanmin(uv(:)) nanmin(vv(:)) nanmin(wv(:))];
err = [];

if (~isempty(vecunits.vel)),
    units = [vecunits.vel{1} '/' vecunits.vel{2}];
else
    units = '';
end;

if ~data.PIV.isStereo
    N = sum(isfinite(uv) & isfinite(vv));
else
    N = sum(isfinite(uv) & isfinite(vv) * isfinite(wv));
end

ctr = [nanmean(xd(:)) nanmean(yd(:))];

% -------------------------------------------------
function [mn,err,units,N,ctr] = apCalcMax(data,rgn,xd,yd,ud,vd, wd, ...
                                   xv,yv,uv,vv,wv, vecunits)

mn = [nanmax(uv(:)) nanmax(vv(:)) nanmax(wv(:))];
err = [];

if (~isempty(vecunits.vel)),
    units = [vecunits.vel{1} '/' vecunits.vel{2}];
else
    units = '';
end;

if ~data.PIV.isStereo
    N = sum(isfinite(uv) & isfinite(vv));
else
    N = sum(isfinite(uv) & isfinite(vv) * isfinite(wv));
end

ctr = [nanmean(xd(:)) nanmean(yd(:))];

% -------------------------------------------------
function [mn,err,units,N,ctr] = apCalcMaxMag(data,rgn,xd,yd,ud,vd,wd, ...
                                   xv,yv,uv,vv,wv, vecunits)
if ~data.PIV.isStereo,
    mn = nanmax(sqrt(uv.^2 + vv.^2));
else
    mn = nanmax(sqrt(uv.^2 + vv.^2 + wv.^2));
end
err = [];

if (~isempty(vecunits.vel)),
    units = [vecunits.vel{1} '/' vecunits.vel{2}];
else
    units = '';
end;

if ~data.PIV.isStereo
    N = sum(isfinite(uv) & isfinite(vv));
else
    N = sum(isfinite(uv) & isfinite(vv) & isfinite(wv));
end

ctr = [nanmean(xd(:)) nanmean(yd(:))];

% -------------------------------------------------
function [flux,err,units,N,ctr] = apCalcFlux(data,rgn,xd,yd,ud,vd,wd, ...
                                   xv,yv,uv,vv,ww, vecunits)

% get the scale factor between position units and velocity length unit
[posscalefac,unitserror] = feval(data.generalFcns.apConvertUnits, ...
                                 vecunits.pos,vecunits.vel{1});

if (~unitserror),
    vecunits.pos = vecunits.vel{1};
end;

xd = xd.*posscalefac;
yd = yd.*posscalefac;

s = [0 cumsum(sqrt(diff(xd).^2 + diff(yd).^2))];

dxds = deriv(s,xd);
dyds = deriv(s,yd);
mag = sqrt(dxds.^2 + dyds.^2);

%normal vector
nx = -dyds./mag;
ny = dxds./mag;

%volume flux
vol = rho * trapz(s,ud.*nx + vd.*ny);
ke = rho * trapz(s,(ud.^2 + vd.^2)*(ud.nx + vd.ny));
momx = rho * trapz(s,ud .* (ud.*nx + vd.*ny));
momy = rho * trapz(s,vd .* (ud.*nx + vd.*ny));

flux.volume = vol;
flux.ke = ke;
flux.momx = momx;
flux.momy = momy;

N = sum(isfinite(uv) & isfinite(vv));
ctr = [nanmean(xd(:)) nanmean(yd(:))];

units = {'','','',''};
if (~isempty(vecunits.pos) && ~isempty(vecunits.vel)),
    if (strcmp(vecunits.pos,vecunits.vel{1})),
        units{1} = sprintf('kg/%s/%s)',vecunits.vel{2},vecunits.vel{1});
        units{2} = sprintf('kg.%s^2/%s^3/%s',vecunits.vel{1},vecunits.vel{2},vecunits.vel{1});
        units{3} = sprintf('kg.%s/%s^2/%s',vecunits.vel{1},vecunits.vel{2},vecunits.vel{1});
        units{4} = units{3};
    else
        warning('Flux units problem');
    end;
end;
err = [];

% -------------------------------------------------
function [x,y,u,v,w, mask] = apGetVecInside(data,xd,yd)

if (size(data.PIV.x,3) > 1),
    x = data.PIV.x(:,:,data.curFrame);
    y = data.PIV.y(:,:,data.curFrame);
else
    x = data.PIV.x;
    y = data.PIV.y;
end;

% fully close the curve
xd(end+1) = xd(1);
yd(end+1) = yd(1);

mask = inpolygon(x,y, xd,yd);

if isfield(data.PIV,'us')
    u = data.PIV.us(:,:,data.curFrame) - data.subtractVector(1);
    v = data.PIV.vs(:,:,data.curFrame) - data.subtractVector(2);
else
    u = data.PIV.u(:,:,data.curFrame) - data.subtractVector(1);
    v = data.PIV.v(:,:,data.curFrame) - data.subtractVector(2);
end
if ~data.PIV.isStereo, 
    w = []; 
elseif isfield(data.PIV,'ws')
    w = data.PIV.ws(:,:,data.curFrame) - data.subtractVector(3);
else
    w = data.PIV.w(:,:,data.curFrame) - data.subtractVector(3);
end
u = u(mask);
v = v(mask);
x = x(mask);
y = y(mask);
if data.PIV.isStereo, 
    w = w(mask); 
end % otherwise w = []


