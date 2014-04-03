function movieplot(varargin)
% function moviePlot(x,y,opts,...)
% Takes more or less the same parameters as plot, but if any of the
% parameters have multiple columns, it plots each column as a separate
% frame.  Note that all multiple column parameters must have the *same*
% number of columns, otherwise it won't recognize them.  Also, the first
% parameter must be multiple column.
%
% If the movie runs too fast, you can pass an option 'delay' followed by
% number of milliseconds to delay after each frame (20-50 is usually reasonable).
%
% Example:
%	size(x) = [4 200];		size(y) = [4 200];
%
% moviePlot(x,y,'k-',x(1,:),y(1,:),'ro',[0; 1],[0 0],'b-','delay',40);
%
% Shows x and y as a movie with 200 frames, plotting a red o at the first
% point in x and y at each frame, and showing a line along the x axis in
% each frame.
%
% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2010, Eric Tytell

opt.delay = 0;
opt.hideaxis = false;
opt.axisequal = true;
opt.axislimits = [];
opt.outputfile = '';
opt.fps = 5;
opt.plotfcn = @plot;
opt.size = [];          %size of AVI image in pixels
opt.background = 'w';

if ((numel(varargin{1}) == 1) && ishandle(varargin{1}) && ...
        strcmp(get(varargin{1},'type'), 'figure'))
    fig = varargin{1};
    p = 2;
else
    fig = gcf;
    p = 1;
end

set(fig,'DoubleBuffer','on', 'Color',opt.background);

[opt,params] = parsevarargin(opt,varargin(p:end),1,'allowno','leaveunknown');

if (length(params) < 2)
    error('Too few arguments');
end;

if (size(params{1},2) == 1),
    N = size(params{2},2);
else
    N = size(params{1},2);
end;

if (~isempty(opt.size))
    set(fig, 'WindowStyle','normal', 'Units','pixels');
    pos = get(fig, 'Position');
    pos([3 4]) = opt.size;
    set(fig, 'Position',pos);
    
    clf;
    if (opt.hideaxis)
        axes('Position',[0 0 1 1]);
    end
end

if (isempty(opt.axislimits)),
    ax = [Inf -Inf Inf -Inf];
    i = 1;
    done = false;
    while (~done && (i <= length(params))),
        if (i+1 > length(params)),
            error('Cannot parse plot parameters');
        end;

        if (ischar(params{i}))
            done = true;
        else        
            ax(1) = min([ax(1); params{i}(:)]);
            ax(2) = max([ax(2); params{i}(:)]);
            ax(3) = min([ax(3); params{i+1}(:)]);
            ax(4) = max([ax(4); params{i+1}(:)]);

            i = i+2;

            if ((i <= length(params)) && ischar(params{i})),
                if (matchlinespec(params{i})),
                    i = i+1;
                end;
            end;
        end
    end;

    opt.axislimits = ax;
end;

if (~isempty(opt.outputfile)),
    aviobj = VideoWriter(opt.outputfile);
    set(aviobj,'FrameRate',opt.fps);
    open(aviobj);
end;

p1 = cell(size(params));
for i = 1:N,
    for j = 1:length(params),
        if (isnumeric(params{j}) && (size(params{j},2) == N)),
            p1{j} = params{j}(:,i);
        else
            p1{j} = params{j};
        end;
    end;
    feval(opt.plotfcn,p1{:});

    if (opt.axisequal),
        axis equal;
    end;
    
    axis(opt.axislimits);
    
    if (opt.hideaxis),
        axis off;
    end;
    
    drawnow;
    
    if (~isempty(opt.outputfile)),
        F = getframe(fig);
        sz = size(F.cdata);
        sz = floor(sz(1:2)/4)*4;
        F.cdata = F.cdata(1:sz(1),1:sz(2),:);
        
        writeVideo(aviobj,F);
    end;
        
    if (opt.delay > 0),
        pause(opt.delay);
    end;
end;

if (~isempty(opt.outputfile)),
    close(aviobj);
end;

