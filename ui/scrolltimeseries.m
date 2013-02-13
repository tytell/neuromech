function scrolltimeseries(varargin)
% function scrolltimeseries(t,y,linspec,...)
% Takes basically the same form of parameters as plot.  Can also take a set of t values
% without a y value, and it will plot vertical lines at those t values.  If it has trouble
% identifying a vertical line t values, pass an empty matrix for y.
% Options: 'totalrange','displayrange','on','off','modal'

i = 1;
plotargs = cell(3,0);
plotopt = {};
vertplotargs = cell(3,0);
opt.totalrange = [-Inf Inf];
opt.displayrange = 50;
opt.off = false;
opt.on = false;
opt.modal = false;
while (i <= length(varargin)),
    if (isnumeric(varargin{i}) && (i < nargin) && isnumeric(varargin{i+1}) && ...
        (ndims(varargin{i}) == 2) && (ndims(varargin{i+1}) == 2) && ...
        all(size(varargin{i}) == size(varargin{i+1}))),

        %two numeric parameters with matching sizes
        a = size(plotargs,2)+1;
        plotargs(1:2,a) = varargin(i:i+1)';
        if ((i+2 <= nargin) && ischar(varargin{i+2}) && matchlinespec(varargin{i+2})),
            plotargs{3,a} = varargin{i+2};
            i = i+3;
        else
            i = i+2;
        end;
    elseif (isnumeric(varargin{i}) && (i < nargin) && isnumeric(varargin{i+1}) && ...
            (ndims(varargin{i}) == 2) && (ndims(varargin{i+1}) == 2) && ...
            any(size(varargin{i}) == 1) && ...
            any(size(varargin{i+1}) == length(varargin{i}))),
        
        %vector time parameter and matrix y parameter
        if (size(varargin{i+1},2) == length(varargin{i})),
            varargin{i+1} = varargin{i+1}';
        end;
        
        a = size(plotargs,2)+1;
        plotargs{1,a} = varargin{i}(:);
        plotargs{2,a} = varargin{i+1};
        
        if ((i+2 <= nargin) && ischar(varargin{i+2}) && matchlinespec(varargin{i+2})),
            plotargs{3,a} = varargin{i+2};
            i = i+3;
        else
            i = i+2;
        end;
    elseif (isnumeric(varargin{i}) && ...
            ((i == nargin) || ...
             ((i < nargin) && isnumeric(varargin{i+1}) && isempty(varargin{i+1})) || ...
             ((i < nargin) && ischar(varargin{i+1})))),
        
        %just a time parameter
        a = size(vertplotargs,2)+1;
        vertplotargs{1,a} = varargin{i};
       
        if ((i+1 <= nargin) && ischar(varargin{i+1}) && matchlinespec(varargin{i+1})),
            vertplotargs{3,a} = varargin{i+1};
            i = i+2;
        end;
    elseif (ischar(varargin{i})),
        switch lower(varargin{i}),
          case {'totalrange','displayrange'},
            opt.(varargin{i}) = varargin{i+1};
            i = i+2;
            
          case {'on','off','modal','toggle'},
            opt.(varargin{i}) = true;
            i = i+1;
            
          case {'linewidth','markeredgecolor','markerfacecolor','markersize',...
                'clipping','color','displayname','tag','userdata'},
            plotopt(end+(1:2)) = varargin(i:i+1);
            i = i+2;
            
          otherwise,
            error('Unrecognized option %s', varargin{i});
        end;
    else
        error('Unrecognized argument at position %d',i);
    end;
end;

if (opt.off),
    set(gcf,'KeyPressFcn',[]);
    return;
elseif (~opt.on && ~isempty(plotargs)),
    for i = 1:size(plotargs,2),
        inrng = (plotargs{1,i} >= opt.totalrange(1)) & ...
                (plotargs{1,i} <= opt.totalrange(2));
        plotargs{1,i} = plotargs{1,i}(inrng);
        if (size(plotargs{1,i},2) == 1),
            plotargs{2,i} = plotargs{2,i}(inrng,:);
        else
            plotargs{2,i} = plotargs{2,i}(inrng);
        end;
    end;
    for i = 1:size(vertplotargs,2),
        inrng = (vertplotargs{1,i} >= opt.totalrange(1)) & ...
                (vertplotargs{1,i} <= opt.totalrange(2));
        vertplotargs{1,i} = vertplotargs{1,i}(inrng);
    end;

    plot(plotargs{:},plotopt{:});
    yl = ylim;
    for i = 1:size(vertplotargs,2),
        tt = vertplotargs{1,i}(:)';
        nt = size(tt,2);
        tt = cat(1,tt([1;1],:),NaN(1,nt));
        vertplotargs{1,i} = tt(:);
        
        yy = [yl(:); NaN];
        yy = yy(:,ones(1,nt));
        vertplotargs{2,i} = yy(:);
    end;

    if (~isempty(vertplotargs)),
        vertplot(vertplotargs{:},plotopt{:});
    end;
end;

if (isfinite(opt.displayrange)),
    xl = xlim;
    xl(2) = xl(1) + opt.displayrange;
    
    xlim(xl);
end;

set(gcf,'KeyPressFcn',@stsKeyPressFcn);

if (opt.modal),
    uiwait;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stsKeyPressFcn(src,event)

xl = xlim;
d = xl(2)-xl(1);

sz = 1;
skip = false;
if (~isempty(event.Modifier)),
    switch event.Modifier{1},
      case 'control',
        sz = 10;
    
      case 'shift',
        sz = 0.1;
     
      case 'alt',
        sz = 1;
        skip = true;
    end;
end;

switch event.Key,
  case {'rightarrow','return',' '},         % right arrow, return or space
    if (skip),
        hln = findobj(gca,'Type','line');
        xd = get(hln,'XData');
        
        isblank = true;
        tnext = zeros(size(xd));
        i = 1;
        while (isblank && (i <= length(xd)))
            tn = first(xd{i}(:),xd{i}(:) >= xl(1));
            if (isempty(tn)),
                tnext(i) = NaN;
            elseif (tn <= xl(2)),
                isblank = false;
            else
                tnext(i) = tn;
            end;
            i = i+1;
        end;
        
        if (isblank && any(~isnan(tnext))),
            xl = min(tnext) + [0 d];
            xlim(xl);
        end;
    else
        xl = xl + 0.9*sz*d;
        xlim(xl);
    end;
    
  case {'leftarrow','backspace'},                        % left arrow
    if (skip),
        hln = findobj(gca,'Type','line');
        xd = get(hln,'XData');
        
        isblank = true;
        tprev = zeros(size(xd));
        i = 1;
        while (isblank && (i <= length(xd)))
            tp = last(xd{i}(:),xd{i}(:) <= xl(end));
            if (isempty(tp))
                tprev(i) = NaN;
            elseif (tp >= xl(1)),
                isblank = false;
            else
                tprev(i) = tp;
            end;
            i = i+1;
        end;
        
        if (isblank && any(~isnan(tprev))),
            xl = max(tprev) + [-d 0];
            xlim(xl);
        end;
    else
        xl = xl - 0.9*sz*d;
        xlim(xl);
    end;
    
  case 'uparrow',
    d = d * (1/(1+sz));
    xl = (xl(1)+xl(2))/2 + [-d d]/2;
    xlim(xl);
    
  case 'downarrow',
    d = d * (1+sz);
    xl = (xl(1)+xl(2))/2 + [-d d]/2;
    xlim(xl);
    
  case 'g',
    ctr = inputdlg('Center on time:','Center');
    xl = ctr + [-d d]/2;
    xlim(xl);
    
  case 'q',
    uiresume;
end;



            