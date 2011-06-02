function varargout = plottimeseries(varargin)
% plottimeseries(plotargs, options)
%
% Plots time series data as multiple rows.  plotargs are any options you
% could pass to the plot function, as long as the first parameter is time.
% Scales the figure so that it fits on a normal page (portrait layout).
% 
% Options (all positions in inches):
%   'margins' Four elment vector with page margins. Format is [left right    
%             top bottom].  (Default [0.5 0.5 0.75 1]
%   'lineheight' Height of the rows. (Default 0.75)
%   'linegap' Gap between rows. (Default 0.2)
%   'linedur' Time duration displayed on each row. (Default 10sec).
%   'ylim' Y limits for each row.
%   'print' Automatically print.
%   'add' Add the plot data to an existing time series plot.

% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2011, Eric Tytell <tytell at jhu dot edu>

margins = [0.5 0.5 0.75 1];         % left right top bottom in inches
lineheight = 0.75;                  % height of a line of traces in inches
linegap = 0.2;                      % gap between lines in inches
linedur = 10;                       % duration of a single line in sec
isPrint = false;
isAdd = false;
usebreaks = false;
plotstyle = 'plot';
ylimdef = [];

if (ischar(varargin{1})),
    switch lower(varargin{1}),
        case {'p','plot'},
            plotstyle = 'plot';
            
        case {'c','colorbar'},
            plotstyle = 'colorbar';
            
        otherwise,
            error('Unrecognized option %s',varargin{1});
    end;
    p = 2;
else     
    p = 1;
end;

a = 1;
t = cell(1,1);
y = cell(1,1);
linestyle = cell(1,1);
plotopts = {};
while (p <= nargin),
    if ((p+1 <= nargin) && isnumeric(varargin{p}) && isnumeric(varargin{p+1})),
        linestyle1 = '';
        step = 2;
        if (~strcmp(plotstyle,'colorbar') && (p+2 <= nargin) && ischar(varargin{p+2})),
            if (regexpi(varargin{p+2}, '^[cmyrgbwk]?(-|--|:|-.)?[.+o*xsd^v><ph]?$')),
                linestyle1 = varargin{p+2};
                step = 3;
            end;
        end;
        if (any(size(varargin{p}) == 1) && any(size(varargin{p+1}) == length(varargin{p}))),            
            t1 = shiftdim(varargin{p});
            if (size(varargin{p+1},1) ~= length(t1)),
                varargin{p+1} = varargin{p+1}';
            end;
            for i = 1:size(varargin{p+1},2),
                t{a} = t1;
                y{a} = varargin{p+1}(:,i);
                linestyle{a} = linestyle1;
                a = a+1;
            end;
            p = p+step;
        elseif ((ndims(varargin{p}) == ndims(varargin{p+1})) && all(size(varargin{p}) == size(varargin{p+1}))),
            if (size(varargin{p},1) < size(varargin{p},2)),
                varargin{p} = varargin{p}';
                varargin{p+1} = varargin{p+1}';
            end;
            
            for i = 1:size(varargin{p+1},2),
                t{a} = varargin{p}(:,i);
                y{a} = varargin{p+1}(:,i);
                linestyle{a} = linestyle1;
                a = a+1;
            end;
            p = p+step;
        else
            error('Unrecognized argument');
        end;
    elseif (ischar(varargin{p})),
        switch lower(varargin{p}),
          case 'margins',
            margins = varargin{p+1};
            p = p+2;
            
          case 'lineheight',
            lineheight = varargin{p+1};
            p = p+2;
            
          case 'linegap',
            linegap = varargin{p+1};
            p = p+2;
            
          case 'linedur',
            linedur = varargin{p+1};
            p = p+2;
            
          case 'print',
            isPrint = true;
            p = p+1;
            
          case 'add',
            isAdd = true;
            p = p+1;
            
          case 'ylim',
            ylimdef = varargin{p+1};
            p = p+2;

          case 'usebreaks',
            usebreaks = true;
            p = p+1;
            
          case 'nobreaks',
            usebreaks = false;
            p = p+1;
            
          otherwise,
            if ((p < nargin) && isnumeric(varargin{p+1})),
                plotopts(end+1:end+2) = varargin(p:p+1);
                p = p+2;
            else
                plotopts{end+1} = varargin{p};
                p = p+1;
            end;
        end;
    else
        error('Unrecognized parameter.');
    end;
end;
        
        
%display the figure on screen at 1.5x smaller than the page printout
figReduce = 1.5;

%create the figure and set up the scaling correctly
fig = findobj('Tag','plottimeseriesfigure');
if (isempty(fig)),
    fig = figure('Tag','plottimeseriesfigure');
    isAdd = false;
elseif (~isAdd),
    close(fig(2:end));
    fig = fig(1);
    figure(fig);
else
    hAxAll = findobj('Tag','plottimeseriesaxis');
    lim = get(hAxAll, 'XLim');
    lim = cat(1,lim{:});
    
    [q,ord] = sort(lim(:,1));
    lim = lim(ord,:);
    hAxAll = hAxAll(ord,:);
    
    parents = get(hAxAll,'Parent');
    parents = cat(1,parents{:});
    
    nAxes = max([sum(parents == min(parents)) sum(parents == max(parents))]);
end;

if (~isAdd),
    clf;

    %set up the figure to have the same proportions as the page size
    set(fig,'PaperUnits','inches','Units','inches');
    papersz = get(fig,'PaperSize');
    paperpos = [0 0 papersz];
    figpos = get(fig,'Position');
    %keep the top of the figure in the same place
    figpos(2) = figpos(2) + figpos(4) - papersz(2)/figReduce;
    %but update it to have a half page size with the right aspect ratio
    figpos(3:4) = papersz/figReduce;

    %figure out how many axes are going to fit
    printsz = papersz - [margins(1)+margins(2) margins(3)+margins(4)];
    nAxes = floor(printsz(2) / (lineheight + linegap));
    linegap = (printsz(2) - (nAxes * lineheight))/(nAxes - 1);

    %create the axes layout
    top = papersz(2) - margins(3);
    axlayout = zeros(nAxes+2, 4);
    for i = 1:nAxes,
        axpos = [margins(1) top - i*lineheight - (i-1)*linegap ...
            printsz(1) lineheight];
        axpos = axpos/figReduce;       
        axlayout(i+1,:) = axpos;
    end;

    %create header and footer axes
    axpos = [margins(1) papersz(2)-margins(3) printsz(1) margins(3)/2];
    axpos = axpos/figReduce;
    axlayout(1,:) = axpos;
    
    axpos = [margins(1) margins(4)/2 printsz(1) margins(4)/2];
    axpos = axpos/figReduce;
    axlayout(end,:) = axpos;

    [fig,hAx,hHead,hFoot] = newtimeseriesfig(fig, figpos,paperpos, axlayout);
    
    if (isPrint),
        pagesetupdlg;
    end;
end;

%get the range of times
tinit = min(cellfun(@min,t));
tfinal = max(cellfun(@max,t));

%check the autoscale for the y axis
if (~isAdd),
    if (isempty(ylimdef)),
        axes(hAx(1));
        ymin = min(cellfun(@min,y));
        ymax = max(cellfun(@max,y));
        plot([0 0],[ymin ymax],'o');
        ylimdef = get(hAx(1),'YLim');
    elseif (ischar(ylimdef)),
        switch lower(ylimdef),
            case 'tight',
                ymin = min(cellfun(@min,y));
                ymax = max(cellfun(@max,y));
                ylimdef = [ymin ymax];
        end;
    end;
else
    if (isempty(ylimdef)),
        ylimdef = get(hAxAll(1,1),'YLim');
    elseif (ischar(ylimdef)),
        switch lower(ylimdef),
            case 'tight',
                ylimdef0 = get(hAxAll(1,1),'YLim');
                
                ymin = min(cellfun(@min,y));
                ymax = max(cellfun(@max,y));
                ylimdef = [min(ymin,ylimdef0(1)) max(ymax,ylimdef0(2))];
        end;
    end;            
end;

%and figure out how many total lines we'll need, then how many pages those
%lines will fit on
tlines = tinit:linedur:tfinal;
tlines(end+1) = Inf;

if (usebreaks),
    n = zeros(length(tlines),length(t));
    for i = 1:length(t),
        n(:,i) = histc(t{i},tlines);
    end;

    goodlines = any(n > 0, 2);
    
    nLines = sum(goodlines);
else
    nLines = length(tlines);
    goodlines = true(length(tlines),1);
end;

tlines = tlines(goodlines);

nPages = ceil(nLines/nAxes);
figs = zeros(nPages,1);

if (isAdd),
    if (nLines > length(hAxAll)),
        error('Cannot overlay plot because time range is different');
    end;
end;

%now run through the pages
P = cell(3,length(t));
axall = 1;
ln = 1;
for page = 1:nPages,
    %now run through all the axes on this page
    for ax = 1:nAxes,
        if (isAdd),
            if (axall > length(hAxAll)),
                continue;
            end;
            curax = hAxAll(axall);
            axall = axall+1;
        else
            curax = hAx(ax);
        end;
        
        axes(curax);
        hold on;
                  
        if (ln < length(tlines)),
            set(curax,'Visible','on');
        else
            %make extra axes invisible
            cla(curax,'reset');
            set(curax,'Visible','off');
            continue;
        end;

        switch plotstyle,
            case 'plot',
                for i = 1:length(t),
                    t1 = t{i};

                    istime = (t1 >= tlines(ln)) & (t1 < tlines(ln+1));
                    P{1,i} = t1(istime);
                    P{2,i} = y{i}(istime);
                    P{3,i} = linestyle{i};
                 end;
                plot(P{:},plotopts{:});
                
            case 'colorbar',
                istime = (t1 >= tlines(ln)) & (t1 < tlines(ln+1));
                
                if (any(istime)),
                    a = [tlines(ln) t1(istime)' tlines(ln+1)];
                    a = a([1;1],:);
                    c = y{1}(istime)';
                    c = c([1;1],[1 1:end end]);
                    b = ylimdef(:);
                    b = b(:,ones(1,size(a,2)));
                    pcolor(a,b,c);
                end;
                shading flat;
        end;

        set(curax, 'Color','none', 'Tag','plottimeseriesaxis', ...
            'YLim',ylimdef, 'XLim',tlines([ln ln+1]));

        ln = ln+1;
    end;

    %reverse the layering order of the axes so that the tick labels on
    %upper axes will be visible
    for ax = nAxes:-1:1,
        axes(curax);
    end;
    
    if (~isAdd),
        axes(hFoot);
        delete(allchild(hFoot));
        text(0.5,0, sprintf('Page %d of %d',page,nPages), ...
            'HorizontalAlignment','center','VerticalAlignment','top');
    end;
    
    if (isPrint),
        printdlg;
    elseif ((page < nPages) && ~isAdd),
        figs(page) = fig;
        [fig,hAx,hHead,hFoot] = newtimeseriesfig(-1, figpos,paperpos, axlayout);
    end;
end;

hAxAll = findobj('Tag','plottimeseriesaxis');
linkaxes(hAxAll,'y');

if (nargout == 1),
    varargout = {hAxAll};
end;

            
function [fig,hAx,hHead,hFoot] = newtimeseriesfig(fig, figpos,paperpos, axlayout)

if (~ishandle(fig)),
    fig = figure('Tag','plottimeseriesfigure');
end;

figure(fig);
clf;

set(fig,'PaperUnits','inches','Units','inches', 'Position',figpos, ...
    'PaperPositionMode','manual','PaperPosition',paperpos, 'Color','w', ...
    'Renderer','painters','InvertHardcopy','off');

hAx = zeros(size(axlayout,1)-2,1);
for i = 2:size(axlayout,1)-1,
    hAx(i-1) = axes('Units','inches', 'Color','none', 'LineWidth',1, ...
        'Tag','plottimeseriesaxis');
    set(hAx(i-1), 'Position',axlayout(i,:), 'Units','normalized');
end;

hHead = axes('Units','inches', 'Color','none', 'LineWidth',1, ...
    'Tag','plottimeseriesheader');
set(hHead, 'Position',axlayout(1,:), 'Units','normalized', 'Visible','off');
hFoot = axes('Units','inches', 'Color','none', 'LineWidth',1, ...
    'Tag','plottimeseriesfooter');
set(hFoot, 'Position',axlayout(end,:), 'Units','normalized', 'Visible','off');
