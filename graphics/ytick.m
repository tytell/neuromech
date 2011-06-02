function ytick(varargin)
% function ytick(tickvals)
%   or     ytick(tickvals, ticklabels)
%   or     ytick(tickval,ticklabel, tickval,ticklabel, ...)
%   or     ytick on
%   or     ytick off
%   or     ytick auto

if ((nargin > 1) && all(ishandle(varargin{1}))),
    ax = varargin{1};
    i = 2;
else
    ax = gca;
    i = 1;
end;

tickopt = {};
tickmodeopt = {};
ticklabelopt = {};

while ((i <= nargin) && ischar(varargin{i})),
    switch lower(varargin{i}),
      case 'on',
        tickmodeopt = {'YTickMode','auto'};
        i = i+1;
        
      case 'off',
        tickopt = {'YTick',[]};
        tickmodeopt = {'YTickMode','manual'};
        i = i+1;
        
      case 'auto',
        tickmodeopt = {'YTickMode','auto'};
        i = i+1;
        
      otherwise,
        error('Unrecognized option %s.', varargin{i});
    end;
end;

if ((i+1 <= nargin) && isnumeric(varargin{i}) && ...
    iscellstr(varargin{i+1})),
    tickopt = {'YTick',varargin{i}};
    ticklabelopt = {'YTickLabel',varargin{i+1}};
    tickmodeopt = {'YTickMode','manual'};
elseif ((i+1 <= nargin) && isnumeric(varargin{i}) && (numel(varargin{i}) == 1) && ...
        ischar(varargin{i+1})),
    if (~all(cellfun(@isnumeric,varargin(i:2:end))) || ...
        ~all(cellfun(@ischar,varargin(i+1:2:end)))),
        error('Bad tick position and label array');
    end;
    
    tickopt = {'YTick',cat(2,varargin{i:2:end})};
    ticklabelopt = {'YTickLabel',varargin(i+1:2:end)};    
    tickmodeopt = {'YTickMode','manual'};
elseif ((i <= nargin) && iscellstr(varargin{i})),
    tickopt = {};
    tickmodeopt = {};
    ticklabelopt = {'XTickLabel',varargin{i}};
elseif ((i <= nargin) && isnumeric(varargin{i})),
    tickopt = {'YTick',varargin{i}};
    tickmodeopt = {'YTickMode','manual'};
end;

set(ax,tickopt{:},tickmodeopt{:},ticklabelopt{:});


           