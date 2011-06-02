function h = circle(varargin)
% function circle(cx,cy, r, linespec, opts...)
%  or  h = circle(ax,cx,cy, r,...)
%
% Plots a circle at cx,cy with radius r.
%
% Options:
%   'nsegments' - Number of segments.  Default = 32.
%
% Also can take any options that plot takes.
%
% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2010, Eric Tytell <tytell at jhu dot edu>

opt.nsegments = 32;
opt.linespec = 'b-';

if ((numel(varargin{1}) == 1) && ishandle(varargin{1}) && ...
    strcmp(get(varargin{1},'Type'),'axes')),
    ax = varargin{1};
    p = 2;
else
    ax = gca;
    p = 1;
end;

[cx,cy,r] = varargin{p:p+2};
p = p+3;

if ((p <= length(varargin)) && matchlinespec(varargin{p}))
    opt.linespec = varargin{p};
    p = p+1;
end;

[opt,rest] = parsevarargin(opt,varargin(p:end), p, 'leaveunknown');

theta = linspace(0,2*pi,opt.nsegments+1);
x = r*cos(theta);
y = r*sin(theta);

h = plot(ax, x+cx, y+cy, opt.linespec, rest{:});
