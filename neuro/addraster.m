function varargout = addraster(varargin)
% function h = addraster(...)
% Adds a raster plot to the existing axes.
% Equivalent to
%   hold on;
%   raster(...)
%   hold off;

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

hold(ax, 'on');
h = raster(varargin{:});

if (~isHeld),
  hold(ax, 'off');
end;

if (nargout == 1)
	varargout{1} = h;
end;
