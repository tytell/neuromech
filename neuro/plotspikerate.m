function varargout = plotspikerate(varargin)
% function plotspikerate(spiket,...)
%   or     plotspikerate(spiket,trialnum,...)
%   or     plotspikerate('rate',t,spikerate,trialnum,...)
%   or     h = plotspikerate(...)
%
% Options:
%    'binsize' - Size of bin for spike rate estimation.  Can be zero, in which case it
%         uses the mean of the forward and reverse spike intervals at each spike.  
%         Default: 0.05
%    'smooth' - Smoothing, in terms of numbers of bins to do a running average over.
%         Default: 3
%    'horizontal' or 'vertical' - Plots the rate bars horizontally (default)
%        or vertically
%    'barwidth' - Width of the rate bars, in fractions of the average
%        between trial numbers.  (default = 0.6).  Also can be 'join', which means
%        to join separate trials
%    'absbarwidth' - Absolute width of the bar.  Overrides 'barwidth'
%    Can also pass most options that fill takes.

%parse numeric parameters
tbin = [];
spikerate = [];
if ((nargin >= 3) && ischar(varargin{1}) && strcmpi(varargin{1},'rate') && ...
    isnumeric(varargin{2}) && isnumeric(varargin{3})),
    [tbin,spikerate] = deal(varargin{2:3});
    if ((nargin >= 4) && isnumeric(varargin{4})),
        y = varargin{4};
        p = 5;
    else
        y = 1:size(tbin,2);
        p = 4;
    end;
elseif ((nargin >= 2) && isnumeric(varargin{1}) && isnumeric(varargin{2}) && ...
    ((size(varargin{1},2) == numel(varargin{2})) || ...
     all(size(varargin{1}) == size(varargin{2})))),
    [spiket,y] = deal(varargin{1:2});
    p = 3;
else
    spiket = varargin{1};
    y = (1:size(spiket,2));
    
    p = 2;
end;

%defaults
ishoriz = true;
barwidth = 0.8;
patchopts = {};
absbarwidth = [];
binsize = 0.05;
smooth = 3;

%check options
while (p <= length(varargin)),
    switch lower(varargin{p}),
      case 'binsize',
        binsize = varargin{p+1};
        p = p+2;
        
      case 'smooth',
        smooth = varargin{p+1};
        p = p+2;
        
      case 'horizontal',
        ishoriz = true;
        p = p+1;
        
      case 'vertical',
        ishoriz = false;
        p = p+1;

      case 'barwidth',
        barwidth = varargin{p+1};
        p = p+2;
        
      case 'absbarwidth',
        absbarwidth = varargin{p+1};
        p = p+2;
        
      case {'annotation', 'clipping', 'edgecolor', 'facecolor', 'displayname', ...
            'linestyle', 'linewidth', 'markeredgecolor','markerfacecolor', ...
            'markersize','tag'},
        % patch options
        patchopts(end+1:end+2) = varargin(p:p+1);
        p = p+2;
        
      otherwise,
        error('Unrecognized option %s', varargin{p});
    end;
end;

%check parameter sizes
if (isempty(spikerate)),
    if (ndims(spiket) > 2),
        error('Cannot handle high dimensional spike data');
    end;

    if (any(size(spiket) ~= size(y))),
        error('y positions (trials) must have the same size as spiket');
    end;
end;

%build the spike rate info
if (isempty(spikerate)),
    [tbin,spikerate] = firingrate(spiket, binsize, smooth);
    iseven = true;
else
    dt = diff(tbin);
    
    if (any(abs(dt(:) - flatten(dt(ones(size(dt,1),1),:))) > 1e-5)),
        iseven = false;
    else
        iseven = true;
    end;
end;

if (size(y,1) == 1),
    dy = nanmedian(diff(sort(y)));
    y = repmat(y,[size(tbin,1) 1]);
else
    dy = nanmedian(flatten(diff(y,[],2)));
end;

if (size(spikerate,2) == 1),           % only one trial
    dy = 1;
end;

if (size(tbin,2) == 1),
    tbin = tbin(:,ones(1,size(spikerate,2)));
end;

%check for absbarwidth
if (isempty(absbarwidth) && ~ischar(barwidth)),
    absbarwidth = barwidth * dy;
elseif (isempty(barwidth)),
    absbarwidth = 0.8 * dy;
end;

if (ishold),
    isheld = true;
    cla;
else
    isheld = false;
    cla reset;
end;

hold on;

if (iseven),
    h = -1 * ones(size(spikerate,2),1);
    for i = 1:size(spikerate,2),
        I = repmat(spikerate(:,i),[1 2]);
        
        good = isfinite(tbin(:,i));
        xx = [first(tbin(:,i),good) last(tbin(:,i),good)];
        if (ischar(barwidth) && strcmpi(barwidth,'join')),
            if ((i > 1) && (i < size(spikerate,2))),
                yy = (y([1 1],i) + y(1,[i-1 i+1])')/2;
            elseif (i == 1),
                yy = [y(1,i)-dy/2 (y(1,i)+y(1,i+1))/2];
            elseif (i == size(spikerate,2)),
                yy = [(y(1,i)+y(1,i-1))/2 y(1,i)+dy/2];
            end;
        else
            yy = y(1,i) + [-1 1]*absbarwidth / 4;
        end;
        
        if (all(isfinite(xx)) && all(isfinite(yy))),
            if (ishoriz),
                h(i) = image('XData',xx','YData',yy','CData',I', ...
                             'CDataMapping','scaled');
            
            else
                h(i) = image('XData',yy,'YData',xx,'CData',I, ...
                             'CDataMapping','scaled');
            end;
        end;
    end;
    
    set(gca,'CLim',[min(spikerate(:)) max(spikerate(:))]);
else
    error('Cannot handle unevenly spaced bins yet.');
end;

if (~isheld),
    hold off;
end;

if (nargout == 1),
    varargout = {h};
end;
