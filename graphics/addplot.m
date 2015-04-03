function varargout = addplot(varargin)
% function h = addplot(...)
% Adds a plot to the existing axes.
% Equivalent to
%   hold on;
%   plot(...)
%   hold off;
% In R2014b and higher, defaults to the older standard of starting the
% color order over again, each time addplot is called.  If you want the
% colors to continuing cycling, use
%   addplot(..., 'restartorder',false);
%
% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2010, Eric Tytell <tytell at jhu dot edu>

opt.restartorder = true;
[opt,args] = parsevarargin(opt,varargin,'leaveunknown','exact');

if ((numel(args{1}) == 1) && ishandle(args{1}) && ...
    strcmp(get(args{1},'Type'),'axes')),
    ax = args{1};
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

if ~verLessThan('matlab','8.4') && opt.restartorder
    % deal with the fact that 
    set(ax,'ColorOrderIndex',1);
end
    
hold(ax, 'on');
h = plot(args{:});

if (~isHeld),
  hold(ax, 'off');
end;

if (nargout == 1)
	varargout{1} = h;
end;
