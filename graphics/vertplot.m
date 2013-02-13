function varargout = vertplot(varargin)
% function vertplot(x,linestyle...)
% Adds vertical lines at particular x values

% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2011, Eric Tytell <tytell at jhu dot edu>

if ((numel(varargin{1}) == 1) && ...
        ishandle(varargin{1}) && strcmp(get(varargin{1},'Type'),'axes')),
    ax = varargin{1};
    i = 2;
else
    ax = gca;
    i = 1;
end;

%check the hold status of the plot
switch get(ax,'NextPlot'),
 case 'replace',
  isHeld = false;
 case 'replacechildren',
  isHeld = false;
 otherwise
  isHeld = true;
end;

yl = makecol(get(ax,'YLim'));

a = 1;
done = false;
opts = {};
while (~done && (i <= length(varargin))),
    if (isnumeric(varargin{i})),
        x1 = varargin{i}(:);
        plotarg{1,a} = repmat(makerow(x1),[2 1]);
        plotarg{2,a} = repmat(yl,[1 length(x1)]);
        
        if ((i < length(varargin)) && matchlinespec(varargin{i+1})),
            plotarg{3,a} = varargin{i+1};
            i = i+1;
        else
            plotarg{3,a} = '';
        end;
        i = i+1;
        a = a+1;
    else
        opts = varargin(i:end);
        done = true;
    end;
end;

hold(ax,'on');

h = plot(ax,plotarg{:},opts{:});

if (~isHeld),
  hold(ax,'off');
end;

if (nargout == 1)
	varargout{1} = h;
end;
