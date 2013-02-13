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

set(gcf,'DoubleBuffer','on');

[opt,p] = parsevarargin(opt,varargin,1,'allowno','leaveunknown');

if (length(p) < 2)
    error('Too few arguments');
end;

if (size(p{1},2) == 1),
    N = size(p{2},2);
else
    N = size(p{1},2);
end;

if (isempty(opt.axislimits)),
    ax = [Inf -Inf Inf -Inf];
    i = 1;
    done = false;
    while (~done && (i <= length(p))),
        if (i+1 > length(p)),
            error('Cannot parse plot parameters');
        end;

        ax(1) = min([ax(1); p{i}(:)]);
        ax(2) = max([ax(2); p{i}(:)]);
        ax(3) = min([ax(3); p{i+1}(:)]);
        ax(4) = max([ax(4); p{i+1}(:)]);
        
        i = i+2;
        
        if ((i <= length(p)) && ischar(p{i})),
            if (matchlinespec(p{i})),
                i = i+1;
            else
                %this should be options for the plot function
                done = true;
            end;
        end;
    end;

    opt.axislimits = ax;
end;

if (~isempty(opt.outputfile)),
    aviobj = avifile(opt.outputfile,'fps',opt.fps);
end;

fig = gcf;

p1 = cell(size(p));
for i = 1:N,
    for j = 1:length(p),
        if (isnumeric(p{j}) && (size(p{j},2) == N)),
            p1{j} = p{j}(:,i);
        else
            p1{j} = p{j};
        end;
    end;
    plot(p1{:});

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
        
        aviobj = addframe(aviobj,F);
    end;
        
    if (opt.delay > 0),
        pause(opt.delay);
    end;
end;

if (~isempty(opt.outputfile)),
    aviobj = close(aviobj);
end;

