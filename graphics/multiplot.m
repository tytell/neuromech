function varargout = multiplot(varargin)
% function multiplot(...)
%
% Overlays multiple types of plots and allows easy stepping through multiple
% frames.  Takes any number of cell arrays that specify the plotting
% function to be called and its usual parameters.
%
% Example:
%    multiplot({'imagesc',ix,iy,I},{'plot',a,b});
%
% The initial element in each cell array is the name or a function handle
% for a graphics function.  Functions also have nicknames, as follows:
%
%   'plot', 'p'
%   'imagesc', 'v'
%   'imshow', 'i'
%   'quiverc', 'q'
%   'pcolor', 'b'
%   'contourf', 'c'
%   show movie, 'm' (takes a movie file name, or x,y, and a file name)
%
% If any of the elements are 3D, multiplot becomes interactive and allows
% stepping through multiple frames using the left and right arrow keys.
%
% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2010, Eric Tytell

opt.dim = 3;

if (~iscell(varargin{1})),
    error('Wrong syntax for multiplot.  Use cell arrays.');
end;

[opt,fcns] = parsevarargin(opt,varargin,'leaveunknown');

if (any(~cellfun(@iscell,fcns))),
    error('multiplot:cellargs','All arguments must be cell arrays');
end;

nfcns = length(fcns);

sz = cell(1,nfcns);
for i = 1:nfcns,
    f = fcns{i};

    sz1 = zeros(3,numel(f));
    sz1(1,:) = cellfun('size',f,1);
    sz1(2,:) = cellfun('size',f,2);
    sz1(3,:) = cellfun('size',f,3);

    iscellarg = cellfun(@iscell,f);
    sz1(3,iscellarg) = sz1(1,iscellarg).*sz1(2,iscellarg);
    
    sz{i} = sz1;
end;
sz = cat(2,sz{:});

rem = [];
post = '';
pre = '';
frnum = [];
for i = 1:nfcns,
    if (ischar(fcns{i}{1})),
        switch lower(fcns{i}{1}),
            case 'post',
                post = fcns{i}{2};
                rem(end+1) = i;
            case 'pre',
                pre = fcns{i}{2};
                rem(end+1) = i;
            case {'frame','fr'},
                frnum = fcns{i}{2};
                rem(end+1) = i;
            case {'quiverc','q'},
                fcns{i}{1} = @quiverc;
            case {'imshow','i'},
                fcns{i}{1} = @mpImShow;
            case {'imagesc','v'},
                fcns{i}{1} = @imagesc;
                if ((length(fcns{i}) == 4) && ...
                        ((isnumeric(fcns{i}{2}) && any(isnan(fcns{i}{2}))) || ...
                        (isnumeric(fcns{i}{3}) && any(isnan(fcns{i}{3}))))),
                    warning('imagesc does not handle NaNs in x and y well.');
                end;
            case {'pcolor','b'},
                fcns{i}{1} = @pcolor;
            case {'contourf','c'},
                fcns{i}(1) = @contourf;
            case {'plot','p'},
                fcns{i}{1} = @plot;
            case {'movie','m'},
                fcns{i}{1} = @mpMovShow;
                
                if ((length(fcns{i}) >= 2) && ischar(fcns{i}{2})),
                    fn = fcns{i}{2};
                    j = 2;
                elseif ((length(fcns{i}) >= 4) && isnumeric(fcns{i}{2}) && ...
                        isnumeric(fcns{i}{3}) && ischar(fcns{i}{4})),
                    fn = fcns{i}{4};
                    j = 4;
                else
                    error('Unrecognized movie arguments.');
                end;
                
                [~,~,ext] = fileparts(fn);
                k = size(sz,2)+1;
                if ismember(ext,{'.avi','.tif','.cine','.mpg'})
                    reader = VideoReader2(fn);
                    sz(3,k) = reader.NumberOfFrames;

                    fcns{i}{j} = reader;
                else
                    error('Unrecognized movie file type');
                end;
                if (length(fcns{i}) >= j+1)
                    frames = fcns{i}{j+1};
                    sz(3,k) = length(frames);
                else
                    frames = 1:reader.NumberOfFrames;
                    fcns{i}{j+1} = frames;
                end
                
            otherwise,
                if ((i == 1) && ischar(fcns{i}{1}) && (length(fcns{i}) == 1)),
                    pre = fcns{i}{1};
                    rem(end+1) = i;
                elseif ((i == nfcns) && ischar(fcns{i}{1}) && ...
                        (length(fcns{i}) == 1)),
                    post = fcns{i}{1};
                    rem(end+1) = i;
                else
                    error('Unrecognized function %s.',fcns{i}{1});
                end;
        end;
    end;
end;

switch get(gca,'NextPlot'),
    case 'replace',
        cla reset;
    case 'replacechildren',
        cla;
end;
zoom off;

keep = setdiff(1:nfcns,rem);
fcns = fcns(keep);
nfcns = length(fcns);

data.fcns = fcns;
data.nfcns = nfcns;
data.pre = pre;
data.post = post;
data.fr = 1;
data.handles = cell(1,nfcns);
data.Figure = gcf;

if (any(sz(opt.dim,:) ~= 1)),
    if (isempty(frnum)),
        % calculate the most frequent number of frames that isn't 1
        nfr = sz(opt.dim,:);
        nfr = mode(nfr(nfr ~= 1));
        nind = nfr;

        frnum = 1:nfr;
        frind = 1:nfr;
    else
        [frnum,q,frind] = unique(frnum);
        nfr = length(frnum);
        nind = length(frind);
    end;

    data.nfr = nfr;
    data.nind = nind;
    data.frnum = frnum;
    data.frind = frind;

    if (opt.dim ~= 3),
        switch opt.dim,
            case 1,
                pmt = [2 3 1];
            case 2,
                pmt = [1 3 2];
        end;
        for f = 1:nfcns,
            for i = 1:length(fcns{f}),
                if (size(fcns{f}{i},opt.dim) == nfr),
                    fcns{f}{i} = permute(fcns{f}{i},pmt);
                end;
            end;
        end;
        data.fcns = fcns;
    end;
    
    guidata(data.Figure, data);

    set(data.Figure, 'KeyPressFcn',@mpKeyPress,...
        'DoubleBuffer','on');

    try
        mpDraw(data.Figure, data, 1);
        uiwait(data.Figure);
    catch
        mpCleanup(data.Figure, data);
        rethrow(lasterror);
    end;
    mpCleanup(data.Figure, data);

    if (nargout == 1),
        if (ishandle(data.Figure)),
            data = guidata(data.Figure);
            varargout{1} = data.fr;
        else
            varargout{1} = NaN;
        end;
    end;
else
    data.nfr = 1;
    data.frind = 1;
    data.nind = NaN;

    data = mpDraw(data.Figure, data, 1);
    if (nargout),
        varargout{1} = cat(1,data.handles{:});
    end;
end;

% ---------------------------------------------------
function data = mpDraw(fig, data, fr)

figure(fig);
cla;

hold on;
if (~isempty(data.pre)),
    eval(data.pre);
end;

pgind = find(data.frind == fr);
data.fr = fr;
guidata(fig, data);

for i = 1:data.nfcns,
    for pgi = 1:length(pgind),
        pg = pgind(pgi);

        f = data.fcns{i};
        for j = 1:length(f),
            if (iscell(f{j})),
                if (length(f{j}) == data.nind),
                    f{j} = f{j}{pg};
                else
                    error(['Cell array length does not match number of ' ...
                        'frames.']);
                end;
            elseif (size(f{j},3) == data.nind),
                f{j} = f{j}(:,:,pg);
            elseif (size(f{j},3) == 1),
                if (size(f{j},2) == data.nind),
                    f{j} = f{j}(:,pg);
                elseif (size(f{j},1) == data.nind),
                    f{j} = f{j}(pg,:);
                end;
            end;
        end;
        data.handles{i} = cat(2,data.handles{i},feval(f{:}));
    end;
end;
if (~isempty(data.post)),
    eval(data.post);
end;
hold off;

drawnow;

set(fig, 'Name', sprintf('multiplot Frame %d/%d',fr,data.nfr));

guidata(fig, data);

% ---------------------------------------------------
function h = mpImShow(varargin)

arg = varargin;
if (ischar(arg{1})),
    I = imread(arg{1});
    arg{1} = I;
elseif ((nargin >= 3) & isnumeric(arg{1}) & isnumeric(arg{2}) & ...
        ischar(arg{3})),
    I = imread(arg{3});
    arg{3} = I;
end;

h = imshow6(arg{:},'n');


% ---------------------------------------------------
function h = mpMovShow(varargin)

arg = varargin;
reader = arg{end-1};
fr1 = arg{end};

I = im2double(read(reader, fr1));

h = imshow6(arg{1:end-2},I,'n');

% ---------------------------------------------------
function mpCleanup(fig, data)

if (ishandle(fig)),
    set(fig, 'KeyPressFcn',[],'UserData',[]);
end;

% ---------------------------------------------------
function mpKeyPress(hObj, eventdata)

data = guidata(hObj);
c = get(data.Figure, 'CurrentCharacter');

draw = 0;
switch lower(c),
    case char(28),                         % left arrow
        if (data.fr > 1),
            data.fr = data.fr-1;
            draw = 1;
        end;
    case {char(29),char(13),' '},          % right arrow, return
        if (data.fr < data.nfr),
            data.fr = data.fr+1;
            draw = 1;
        end;
    case 'g',                              % go to frame
        fr = inputdlg('Frame?','Go to frame');
        fr = str2num(fr{1});
        if ((fr >= 1) & (fr <= data.nfr)),
            data.fr = fr;
            draw = 1;
        end;
    case 'q',                              % quit
        uiresume;
end;

if (draw),
    mpDraw(data.Figure, data, data.fr);
end;




