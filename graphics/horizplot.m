function varargout = horizplot(varargin)
% function horizplot(y,linestyle...)
% Adds horizontal lines at particular y values
%
% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2010, Eric Tytell <tytell at jhu dot edu>

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

%get the x limits
xl = makecol(get(ax,'XLim'));

a = 1;
done = false;
opts = {};
while (~done && (i <= length(varargin))),
    if (isnumeric(varargin{i})),
        y1 = varargin{i}(:);
        plotarg{1,a} = repmat(xl,[1 length(y1)]);
        plotarg{2,a} = repmat(makerow(y1),[2 1]);

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
