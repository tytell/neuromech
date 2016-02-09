function [x,y] = manualTrackLine(filename,varargin)

opt.frameskip = 1;
opt.npts = [];
opt.oldpts = {};
opt.trackmode = 'pointfirst';

opt.timetrace = true;
opt.ntimetrace = 5;

opt = parsevarargin(opt,varargin,2);

if ~opt.timetrace
    opt.ntimetrace = 0;
end

if ~isempty(opt.oldpts)
    oldx = opt.oldpts{1};
    oldy = opt.oldpts{2};
else
    oldx = [];
    oldy = [];
end

if (isempty(filename)),
    [fn,pn] = uigetfile('*.avi','Open movie file');
    filename = fullfile(pn,fn);
end;

if iscell(filename)
    reader = [];
    filenames = filename;
    nfr = length(filenames);
else
    reader = VideoReader2(filename);
    filenames = [];
    nfr = reader.NumberOfFrames;
end

curfr = 1;
curpt = 1;

clf;

if ~isempty(filenames)
    I = imread(filenames{curfr});
else
    I = read(reader,curfr);
end
hImage = imshow(I,'InitialMagnification','fit');

data.hPrev(1) = line('XData',[], 'YData',[], 'Marker','x', ...
    'Color','r', 'MarkerSize',12, 'ButtonDownFcn', {@mtlLineButtonDown,1});
data.hPrev(2) = line('XData',[], 'YData',[], 'Marker','o', ...
    'Color','r', 'ButtonDownFcn', {@mtlLineButtonDown,2});
if (opt.ntimetrace > 0)
    switch opt.trackmode,
        case 'framefirst',
            for j = 1:opt.ntimetrace,
                data.hPrev(j+2) = line('XData',[],'YData',[], 'Marker','.', ...
                    'Color','y', 'ButtonDownFcn', {@mtlLineButtonDown,j+2});
            end;
        case 'pointfirst',
            data.hPrev(3) = line('XData',[], 'YData',[], 'Marker','.', ...
                'Color','y', 'ButtonDownFcn', {@mtlLineButtonDown,3});
    end;
end;

data.hTitle = title('');

data.Figure = gcf;
data.Axes = gca;
data.hImage = hImage;
data.reader = reader;
data.filename = filename;
data.filenames = filenames;
data.frameskip = opt.frameskip;
data.trackmode = opt.trackmode;
data.ntimetrace = opt.ntimetrace;
data.nfr = nfr;
data.curfr = curfr;
if (~isempty(opt.oldpts)),
    data.rawx = oldx(1:curpt-1, curfr);
    data.rawy = oldy(1:curpt-1, curfr);
else
    data.rawx = [];
    data.rawy = [];
end;
data.curpt = curpt;
data.closing = false;

switch opt.trackmode,
    case 'framefirst',
        data.rawx = NaN(1,nfr);
        data.rawy = NaN(1,nfr);
        
        if (~isempty(oldx)),
            if (size(oldx,2) > nfr),
                oldx = oldx(:,1:nfr);
                oldy = oldy(:,1:nfr);
            end;
            for fr = 1:size(oldx,2),
                good = isfinite(oldx(:,fr)) & isfinite(oldy(:,fr));
                data.rawx{fr} = oldx(good,fr);
                data.rawy{fr} = oldy(good,fr);
            end;
        end;
        
    case 'pointfirst',
        data.rawx = NaN(1,nfr);
        data.rawy = NaN(1,nfr);
        
        if (~isempty(oldx)),
            if (size(oldx,2) <= nfr),
                data.rawx(1:size(oldx,1),1:size(oldx,2)) = oldx;
                data.rawy(1:size(oldx,1),1:size(oldx,2)) = oldy;
            else
                oldx = oldx(:,1:nfr);
                oldy = oldy(:,1:nfr);
                data.rawx(1:size(oldx,1),:) = oldx;
                data.rawy(1:size(oldy,1),:) = oldy;
            end;
        end;
end;


mtlUpdate(data,{'image','points'});

set(data.Figure, 'DoubleBuffer','on', 'KeyPressFcn', @mtlKeyPress, ...
    'CloseRequestFcn', @mtlCloseFigure);
set(data.hImage, 'ButtonDownFcn', @mtlImageButtonDown);

guidata(data.Figure,data);
uiwait(data.Figure);

data = guidata(data.Figure);
if (data.closing),
    delete(data.Figure);
else
    set(data.Figure, 'KeyPressFcn',[], 'CloseRequestFcn', 'closereq');
    set(data.hImage, 'ButtonDownFcn', []);
    set(data.hPrev, 'ButtonDownFcn', []);
end;

x = data.rawx;
y = data.rawy;

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mtlUpdate(data, what)

if (~iscell(what))
    what = {what};
end;

for i = 1:length(what),
    switch what{i},
        case 'image',
            if ~isempty(data.filenames)
                I = imread(data.filenames{data.curfr});
            else
                I = read(data.reader, data.curfr);
            end
            set(data.hImage, 'CData', I);
            
        case 'points',
            prevframes = data.curfr - (1:data.ntimetrace)*data.frameskip;
            prevframes = prevframes(prevframes >= 1);
            
            set(data.hPrev(2), 'XData', data.rawx(1:data.curpt-1, data.curfr), ...
                'YData', data.rawy(1:data.curpt-1, data.curfr));
                    
            switch data.trackmode,
                case 'framefirst',
                    %nothing here...
                    
                case 'pointfirst',
                    set(data.hPrev(1), 'XData',data.rawx(data.curpt,data.curfr), ...
                        'YData',data.rawy(data.curpt,data.curfr));
                    set(data.hPrev(3), 'XData', data.rawx(data.curpt,prevframes), ...
                        'YData', data.rawy(data.curpt,prevframes));
            end;
    end;
end;

str = sprintf('Frame %d/%d, point %d', data.curfr,data.nfr,data.curpt);
set(data.hTitle,'String',str);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mtlImageButtonDown(obj, eventdata)

data = guidata(obj);

pt = get(data.Axes,'CurrentPoint');

switch (data.trackmode),
    case 'framefirst',
        if (data.curpt > size(data.rawx,1)),
            data.rawx(end+1:data.curpt,:) = NaN;
            data.rawy(end+1:data.curpt,:) = NaN;
        end;
        
        data.rawx(data.curpt,data.curfr) = pt(1,1);
        data.rawy(data.curpt,data.curfr) = pt(1,2);

        data.curpt = data.curpt + 1;
        mtlUpdate(data, 'points');
        
    case 'pointfirst',
        data.rawx(data.curpt,data.curfr) = pt(1,1);
        data.rawy(data.curpt,data.curfr) = pt(1,2);
        
        if (data.curfr <= data.nfr - data.frameskip)
            data.curfr = data.curfr + data.frameskip;
        end;
        mtlUpdate(data, {'points','image'});
end;

guidata(obj,data);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mtlLineButtonDown(obj, eventdata, linenum)

data = guidata(obj);
switch get(data.Figure, 'SelectionType'),
    case 'normal',
        mtlImageButtonDown(obj, eventdata);
        
    case 'alternate',
        % eventually edit the point here, but not yet
end;


function mtlKeyPress(obj,eventdata)

data = guidata(obj);

c = get(obj, 'CurrentCharacter');

switch lower(c),
    case char(13),          % return
        switch (data.trackmode),
            case 'framefirst',
                if (data.curfr <= data.nfr + data.frameskip)
                    data.curfr = data.curfr + data.frameskip;
                end;
                data.curpt = 1;
                
                mtlUpdate(data, {'points','image'});
            case 'pointfirst',
                data.curpt = data.curpt + 1;
                data.curfr = 1;
                if (size(data.rawx,1) < data.curpt),
                    data.rawx(end+1:data.curpt, :) = NaN;
                    data.rawy(end+1:data.curpt, :) = NaN;
                end;
        end;
        mtlUpdate(data, {'points','image'});
                
    case 'n',
        if (strcmp(data.trackmode,'pointfirst'))
            data.curpt = data.curpt + 1;
            if (size(data.rawx,1) < data.curpt),
                data.rawx(end+1:data.curpt, :) = NaN;
                data.rawy(end+1:data.curpt, :) = NaN;
            end;
            mtlUpdate(data, 'points');
        end;
        
    case 'p',
        if (strcmp(data.trackmode,'pointfirst') && (data.curpt > 1)),
            data.curpt = data.curpt - 1;
            mtlUpdate(data, 'points');
        end;
        
    case char(28),              % left arrow
        if (data.curfr > data.frameskip),
            data.curfr = data.curfr - data.frameskip;
            if (strcmp(data.trackmode, 'framefirst')),
                data.curpt = 1;
            end;
               
            mtlUpdate(data, {'image','points'});
        end;
        
    case char(29),              % left arrow
        if (data.curfr <= data.nfr - data.frameskip),
            data.curfr = data.curfr + data.frameskip;
            if (strcmp(data.trackmode, 'framefirst')),
                data.curpt = 1;
            end;
            mtlUpdate(data, {'image','points'});
        end;
        
    case 'q',
        uiresume;
end;

guidata(obj,data);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mtlCloseFigure(fig, eventdata)

data = guidata(fig);

data.closing = true;
guidata(fig,data);

uiresume;

            

