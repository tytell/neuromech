function [x,y] = manualTrackLine(varargin)

frameskip = 1;
npts = [];
oldx = [];
oldy = [];
trackmode = 'pointfirst';
aviname = '';
ntimetrace = [];
ntimetracefr = 1;
ntimetracept = 10;

p = 1;
while (p <= length(varargin)),
    if (~ischar(varargin{p})),
        error('Unrecognized parameter #%d',p);
    end;
    
    switch lower(varargin{p}),
        case 'skip',
            frameskip = varargin{p+1};
            p = p+2;
            
        case 'npts',
            npts = varargin{p+1};
            p = p+2;
            
        case 'oldpts',
            oldx = varargin{p+1};
            oldy = varargin{p+2};
            p = p+3;
            
        case {'framefirst','pointfirst'},
            trackmode = varargin{p};
            p = p+1;

        case 'notimetrace',
            ntimetrace = 0;
            p = p+1;
            
        case 'ntrace',
            ntimetrace = varargin{p+1};
            p = p+2;
            
        otherwise,
            if (exist(varargin{p},'file'))
                aviname = varargin{p};
                p = p+1;
            else
                error('Unrecognized option %s',varargin{p});
            end;
    end;
end;

if (isempty(aviname)),
    [fn,pn] = uigetfile('*.avi','Open movie file');
    aviname = fullfile(pn,fn);
end;

reader = mmreader(aviname);

nfr = reader.NumberOfFrames;

curfr = 1;
curpt = 1;

clf;

I = read(reader,curfr);
hImage = imshow(I,'InitialMagnification','fit');

data.hPrev(1) = line('XData',[], 'YData',[], 'Marker','x', ...
    'Color','r', 'MarkerSize',12, 'ButtonDownFcn', {@mtlLineButtonDown,1});
data.hPrev(2) = line('XData',[], 'YData',[], 'Marker','o', ...
    'Color','r', 'ButtonDownFcn', {@mtlLineButtonDown,2});
if (isempty(ntimetrace)),
    switch trackmode,
        case 'framefirst',
            ntimetrace = ntimetracefr;
            for j = 1:ntimetrace,
                data.hPrev(j+2) = line('XData',[],'YData',[], 'Marker','.', ...
                    'Color','y', 'ButtonDownFcn', {@mtlLineButtonDown,j+2});
            end;
        case 'pointfirst',
            ntimetrace = ntimetracept;
            data.hPrev(3) = line('XData',[], 'YData',[], 'Marker','.', ...
                'Color','y', 'ButtonDownFcn', {@mtlLineButtonDown,3});
    end;
end;

data.hTitle = title('');

data.Figure = gcf;
data.Axes = gca;
data.hImage = hImage;
data.reader = reader;
data.aviname = aviname;
data.frameskip = frameskip;
data.trackmode = trackmode;
data.ntimetrace = ntimetrace;
data.nfr = nfr;
data.curfr = curfr;
if (~isempty(oldx)),
    data.rawx = oldx(1:curpt-1, curfr);
    data.rawy = oldy(1:curpt-1, curfr);
else
    data.rawx = [];
    data.rawy = [];
end;
data.curpt = curpt;
data.closing = false;

switch trackmode,
    case 'framefirst',
        data.rawx = nans(1,nfr);
        data.rawy = nans(1,nfr);
        
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
        data.rawx = nans(1,nfr);
        data.rawy = nans(1,nfr);
        
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
            I = read(data.reader, data.curfr);
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

            

