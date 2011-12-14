function varargout = overlayplot(varargin)
% function h = overlayplot()

if ((nargin >= 1) && (numel(varargin{1}) == 1) && ishandle(varargin{1}) && ...
        strcmp(get(varargin{1},'Type'),'axes')),
    ax = varargin{1};
    p = 2;
else
    ax = gca;
    p = 1;
end;

islinkx = true;
isnew = true;
resetsizes = false;
matchx = [];
matchy = [];
top = [];
axiscolor = [];
while ((p <= nargin) && ischar(varargin{p})),
    switch lower(varargin{p}),
        case {'remove','off'},
            [bot,top] = getoverlayaxes(ax,false);
            if (ishandle(top)),
                delete(top);
            else
                warning('No overlay plot found.');
            end;
            return;
            
        case 'linkx',
            islinkx = true;
            p = p+1;
            
        case 'linky',
            islinkx = false;
            p = p+1;
            
        case 'top',
            [bot,top,isnew] = getoverlayaxes(ax);
            axes(top);
            istop = true;
            p = p+1;
            
        case {'bot','bottom'},
            [bot,top,isnew] = getoverlayaxes(ax);
            axes(bot);
            istop = false;
            p = p+1;
            
        case 'get',
            p = p+1;
            
        case 'matchx',
            if ((p < nargin) && isnumeric(varargin{p+1}) && (numel(varargin{p+1}) == 1)),
                matchx = varargin{p+1};
                p = p+2;
            else
                matchx = 0;
                p = p+1;
            end;

        case 'matchy',
            if ((p < nargin) && isnumeric(varargin{p+1}) && (numel(varargin{p+1}) == 1)),
                matchy = varargin{p+1};
                p = p+2;
            else
                matchy = 0;
                p = p+1;
            end;

        case 'axiscolor',
            axiscolor = varargin{p+1};
            p = p+2;
            
        case 'resetsizes',
            [bot,top,isnew] = getoverlayaxes(ax);
            resetsizes = true;
            p = p+2;
            
        otherwise,
            error('Unrecognized option %s',varargin{p});
    end;
end;

if (isempty(top)),
    [bot,top,isnew,istop] = getoverlayaxes(ax);
end;

if (isnew),
    set([bot top],'Box','off');
    if (islinkx),
        %first turn off the linking
        linkaxes([bot top],'off');
        %and set y to scale automatically
        ylim auto;
        %now clear ticks along the x axis
        set(top,'XTick',[],'YAxisLocation','right');
        %and link the x axes
        linkaxes([bot top],'x');
    else
        linkaxes([bot top],'off');
        set(top,'XTickMode','auto','YTick',[],'XAxisLocation','top');
        linkaxes([bot top],'y');
    end;
    axes(top);
    istop = true;
end;

if (p <= nargin),
    switch get(ax,'NextPlot'),
        case 'replace',
            isHeld = false;
        case 'replacechildren',
            isHeld = false;
        otherwise
            isHeld = true;
    end;
    
    hold on;
    plot(varargin{p:end});
    
    if (~isHeld),
        hold off;
    end;
    
    isnew = false;
end;

if (~isnew),
    if (resetsizes),
        toppos = get(top,'Position');
        botpos = get(bot,'Position');
        set(top,'Position',botpos);
    elseif (~isempty(matchy)),
        botylim = get(bot,'YLim');
        topylim = get(top,'YLim');
        
        %set the top y limits such that matchy shows up at the same
        %position in both graphs
        
        %fractional position of matchy
        fracpos = (matchy - botylim(1))/diff(botylim);
        if ((fracpos < 0) || (fracpos > 1)),
            error('Matching y position %g is not visible in bottom graph',matchy);
        end;
        topfracpos = (matchy - topylim(1))/diff(topylim);
        
        shift = fracpos - topfracpos;
        
        newylim = topylim + shift * diff(topylim);
        set(top,'YLim',newylim);
    elseif (~isempty(matchx)),
        botxlim = get(bot,'XLim');
        topxlim = get(top,'XLim');
        
        %set the top x limits such that matchx shows up at the same
        %position in both graphs
        
        %fractional position of matchx
        fracpos = (matchx - botxlim(1))/diff(botxlim);
        if ((fracpos < 0) || (fracpos > 1)),
            error('Matching x position %g is not visible in bottom graph',matchx);
        end;
        topfracpos = (matchx - topxlim(1))/diff(topxlim);
        
        shift = fracpos - topfracpos;
        
        newxlim = topxlim - shift * diff(topxlim);
        set(top,'XLim',newxlim);
    end;
end;

if (~isempty(axiscolor)),
    if (istop),
        hax = top;
    else
        hax = bot;
    end;
    if (islinkx),
        opt = 'YColor';
    else
        opt = 'XColor';
    end;        
    set(hax,opt,axiscolor);
end;
        
if (nargout == 1)
    if (istop),
        varargout = {top};
    else
        varargout = {bot};
    end;
elseif (nargout == 2),
    varargout = {bot,top};
end;
    
%--------------------------------------------------------------------------
function [b,t,isnew,istop] = getoverlayaxes(h,iscreate)
%find the set of two overlay axes corresponding to the one passed in, which
%might be either the top or bottom axis

isnew = false;
if (nargin == 1),
    iscreate = true;
end;

%get the parent figure
par = get(h,'Parent');

%overlay axes by definition have the same position, so look for the one
%that has the same position as h
pos = get(h,'Position');
ax = get(par,'Children');
h2 = [];
for i = 1:length(ax),
    if ((ax(i) ~= h) && all(get(ax(i),'Position') == pos)),
        if (~isempty(h2)),
            error('Too many overlay axes.');
        else
            h2 = ax(i);
        end;
    end;
end;

if (isempty(h2)),
    if (iscreate),
        %we didn't find anything, so create a top axes for h
        b = h;
        t = axes('Position',pos,'Parent',par);
        
        %set it transparent
        set(t,'Color','none');
        isnew = true;
        istop = true;
    else
        b = h;
        t = -1;
        istop = false;
        isnew = false;
    end;
else
    %we found something, so now figure out which is transparent
    if (strcmp(get(h2,'Color'),'none')),
        b = h;
        t = h2;
        istop = false;
    else
        b = h2;
        t = h;
        istop = true;
    end;
end;
