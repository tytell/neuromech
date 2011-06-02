function xtick(varargin)
% function xtick(tickvals)
%   or     xtick(tickvals, ticklabels)
%   or     xtick(tickval,ticklabel, tickval,ticklabel, ...)
%   or     xtick on
%   or     xtick off
%   or     xtick auto

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
        tickmodeopt = {'XTickMode','auto'};
        i = i+1;
        
      case 'off',
        tickopt = {'XTick',[]};
        tickmodeopt = {'XTickMode','manual'};
        i = i+1;
        
      case 'auto',
        tickmodeopt = {'XTickMode','auto'};
        i = i+1;
        
      otherwise,
        error('Unrecognized option %s.', varargin{i});
    end;
end;

if ((i+1 <= nargin) && isnumeric(varargin{i}) && ...
    iscellstr(varargin{i+1})),
    tickopt = {'XTick',varargin{i}};
    ticklabelopt = {'XTickLabel',varargin{i+1}};
    tickmodeopt = {'XTickMode','manual'};
elseif ((i+1 <= nargin) && isnumeric(varargin{i}) && (numel(varargin{i}) == 1) && ...
        ischar(varargin{i+1})),
    if (~all(cellfun(@isnumeric,varargin(i:2:end))) || ...
        ~all(cellfun(@ischar,varargin(i+1:2:end)))),
        error('Bad tick position and label array');
    end;
    
    tickopt = {'XTick',cat(2,varargin{i:2:end})};
    ticklabelopt = {'XTickLabel',varargin(i+1:2:end)};    
    tickmodeopt = {'XTickMode','manual'};
elseif ((i <= nargin) && iscellstr(varargin{i})),
    tickopt = {};
    tickmodeopt = {};
    ticklabelopt = {'XTickLabel',varargin{i}};
elseif ((i <= nargin) && isnumeric(varargin{i})),
    tickopt = {'XTick',varargin{i}};
    tickmodeopt = {'XTickMode','manual'};
end;

set(ax,tickopt{:},tickmodeopt{:},ticklabelopt{:});


           