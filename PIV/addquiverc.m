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
  isHeld = false;
 case 'replacechildren',
  isHeld = false;
 otherwise
  isHeld = true;
end;

hold(ax, 'on');
hh = quiverc(varargin{:});
if ~isHeld
    hold(ax,'off');
end

if (nargout == 1),
	varargout{1} = hh;
end;

