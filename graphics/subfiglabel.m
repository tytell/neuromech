function varargout = subfiglabel(x,y, label, varargin)

opt.defaultlocation = [0 1];
opt.FontSize = 12;
opt.FontWeight = 'bold';
opt.HorizontalAlignment = 'left';
opt.VerticalAlignment = 'bottom';

if (nargin == 1)
    args = {x};
elseif (nargin == 2)
    args = {x, y};
elseif (nargin >= 3)
    args = [x, y, label, varargin];
end

if ((nargin >= 2) && (numel(args{1}) == 1) && ishandle(args{1}) && ...
        strcmp(get(args{1},'Type'), 'axes'))
    ax = args{1};
    p = 2;
else
    ax = gca;
    p = 1;
end;

if ischar(args{p})
    x = [];
    y = [];
    label = args{p};
    p = p+1;
else
    [x, y, label] = args{p:p+2};
    p = p+3;
end

opt = parsevarargin(opt,args(p:end));

if (isempty(x))
    x = opt.defaultlocation(1);
end
if (isempty(y))
    y = opt.defaultlocation(2);
end

htxt = text(x,y, label, 'Parent',ax, 'FontSize',opt.FontSize, 'FontWeight',opt.FontWeight, ...
    'Units','normalized', ...
    'HorizontalAlignment',opt.HorizontalAlignment, 'VerticalAlignment',opt.VerticalAlignment);

if (nargout == 1)
    varargout = {htxt};
end
