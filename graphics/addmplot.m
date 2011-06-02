function varargout = addmplot(varargin)
% function h = addmplot(...)
% Exactly the same as mplot, but adds the plot to the existing axes.
%
% Mercurial revision hash: $Revision: 18f43cd9074e $ $Date: 2010/08/10 21:11:58 $
% Copyright (c) 2010, Eric Tytell <tytell at jhu dot edu>

if ((numel(varargin{1}) == 1) && ishandle(varargin{1}) && ...
    strcmp(get(varargin{1},'Type'),'axes')),
    ax = varargin{1};
else
    ax = gca;
end;

switch get(ax,'NextPlot'),
 case 'replace',
  isHeld = 0;
 case 'replacechildren',
  isHeld = 0;
 otherwise
  isHeld = 1;
end;

hold(ax,'on');
h = mplot(varargin{:});

if (~isHeld),
  hold(ax,'off');
end;

if (nargout == 1)
	varargout{1} = h;
end;
