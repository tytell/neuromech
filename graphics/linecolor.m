function h = linecolor(x,y,w,c,varargin)
% LINECOLOR  Plots a line with varying color and width
%   h  = linecolor(x,y,w,c,...)
%    or  linecolor(ax,x,y,w,c,...)
%
%   Similar to SCATTER, x and y are the coordinates of the line, w is its
%   width (in image coordinates - NOT points), and c is its color (scaled
%   into the colormap).
%
%   Can take any options that are valid for patch.


if ((numel(x) == 1) && ishandle(x) && strcmp(get(x,'Type'),'axes'))
    ax = x;
    [x,y,w] = deal(y,w,c);
    if ((nargin >= 5) && isnumeric(varargin{1}))
        c = varargin{1};
        varargin = varargin(2:end);
    end;
else
    ax = gca;
end;

switch get(ax,'NextPlot'),
 case 'replace',
  cla(ax,'reset');
 case 'replacechildren',
  cla(ax);
end;

if (size(y,1) == 1)
    y = y';
    x = x';
    w = w';
    c = c';
end;
if (any(size(x) == 1) && (length(x) == size(y,1)))
    x = x(:);
    x = x(:,ones(1,size(y,2)));
end;

if (isempty(w) || (all(w(:) == 0)))
    h = patch(x,y,c,'Parent',ax);
    set(h,'FaceColor','none','EdgeColor','flat', varargin{:});
else
    dx = diff(x);
    dy = diff(y);
    ds = sqrt(dx.^2 + dy.^2);
    
    dxds = dx ./ ds;
    dyds = dy ./ ds;
    dxds = dxds([1:end end],:);
    dyds = dyds([1:end end],:);
    
    x1 = shiftdim(x + w.*dyds, -1);
    x2 = shiftdim(x - w.*dyds, -1);
    y1 = shiftdim(y - w.*dxds, -1);
    y2 = shiftdim(y + w.*dxds, -1);
    
    xx = cat(1,x1(1,1:end-1,:),x2(1,1:end-1,:),x2(1,2:end,:),x1(1,2:end,:));
    yy = cat(1,y1(1,1:end-1,:),y2(1,1:end-1,:),y2(1,2:end,:),y1(1,2:end,:));
    
    c = shiftdim(c,-1);
    c = cat(1,c([1 1],1:end-1,:),c([1 1],2:end,:));
    
    h = zeros(size(xx,3),1);
    for i = 1:size(xx,3)
        h(i) = patch(xx(:,:,i),yy(:,:,i),c(:,:,i),'Parent',ax);
    end;
    set(h,'FaceColor','interp', 'EdgeColor','none', varargin{:});
end;




