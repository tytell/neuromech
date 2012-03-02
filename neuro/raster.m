function varargout = raster(varargin)
% function raster(spiket,...)
%   or     raster(spiket,trialnum,...)
%   or     h = raster(...)
%
% Options:
%    'horizontal' or 'vertical' - Plots the rasters horizontally (default)
%        or vertically
%    'barlength' - Length of the raster bars, in fractions of the average
%        between trial numbers.  (default = 0.6)
%    'absbarlength' - Absolute length of the bar.  Overrides 'barlength'
%    Can also pass most options that plot takes.

%parse numeric parameters
if ((nargin >= 2) && (numel(varargin{1}) == 1) && ishandle(varargin{1})),
    ax = varargin{1};
    p = 2;
else
    ax = gca;
    p = 1;
end;
if ((nargin >= p+1) && isnumeric(varargin{p}) && isnumeric(varargin{p+1}) && ...
    ((size(varargin{p},2) == numel(varargin{p+1})) || ...
     all(size(varargin{p}) == size(varargin{p+1})))),
    [spiket,y] = deal(varargin{p:p+1});
    p = p+2;
else
    spiket = varargin{p};
    y = (1:size(spiket,2));
    
    p = p+1;
end;

%defaults
opt.ishoriz = true;
opt.barlength = 0.6;
opt.absbarlength = [];

[opt,plotopts] = parsevarargin(opt,varargin(p:end),'leaveunknown');

%check parameter sizes
if (size(spiket,2) > size(spiket,1)),
    warning(['Trials appear to be in rows.  They should be columns.  Plot ' ...
             'may look weird.']);
end;
if (ndims(spiket) > 2),
    error('Cannot handle high dimensional spike data');
end;

if (size(y,1) == 1),
    dy = nanmedian(diff(y));
    y = repmat(y,[size(spiket,1) 1]);
else
    uy = unique(y(isfinite(y)));
    dy = range(y(isfinite(y))) / length(uy);
end;

if (size(spiket,2) == 1),           % only one trial
    dy = 1;
end;

if (any(size(spiket) ~= size(y))),
    error('y positions (trials) must have the same size as spiket');
end;

%check for absbarlength
if (isempty(opt.absbarlength)),
    absbarlength = opt.barlength * dy;
else
    absbarlength = opt.absbarlength;
end;

if (~isempty(plotopts) && ischar(plotopts{1}))
    [islinespec,~,mark,linestyle] = matchlinespec(plotopts{1});
else
    islinespec = false;
end;
if (islinespec && ~isempty(mark) && isempty(linestyle))
    yy = y(:);
    xx = spiket(:);
else
    %build up the x and y coordinates
    yy = repmat(shiftdim(y,-1),[3 1]);
    yy = yy + repmat([-1; 1; NaN]*absbarlength/2, [1 size(yy,2) size(yy,3)]);
    yy = reshape(yy,[3*size(yy,2) size(yy,3)]);
    
    xx = repmat(shiftdim(spiket,-1), [3 1]);
    xx = reshape(xx,[3*size(xx,2) size(xx,3)]);
end;

%plot horizontally or vertically
if (opt.ishoriz),
    h = plot(ax, xx,yy, plotopts{:});
else
    h = plot(ax, yy,xx, plotopts{:});
end;

if (nargout == 1),
    varargout = {h};
end;
