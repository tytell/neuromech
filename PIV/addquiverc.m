function varargout = addquiverc(varargin)

charargs = find(cellfun(@ischar,varargin));
axopt = find(strcmpi(varargin(charargs),'axes'));
axopt = charargs(axopt); 
if (~isempty(axopt) && (axopt < nargin)),
    ax = varargin{axopt+1};
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
hh = quiverc(varargin{:});

if (nargout == 1),
	varargout{1} = hh;
end;

