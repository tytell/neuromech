function htxt = labellines(labels, varargin)

opts.location = 'best';
opts.rotation = 0;

hax = [];
hln = [];
if (all(ishandle(labels)) && (nargin >= 2))
    h = labels;
    if ((numel(h) == 1) && strcmp(get(h,'Type'),'axes'))
        hax = h;
    else
        hln = h;
    end
    labels = varargin{1};
    p = 2;
else
    p = 1;
end

opts = parsevarargin(opts,varargin(p:end), p+1, 'typecheck',false);

if isempty(hax)
    if ~isempty(hln)
        hax = get(hln(1),'Parent');
    else
        hax = gca;
    end
end
if isempty(hln)
    hln = findobj(hax, 'Type','line');
    %findobj always gives us the object in the reverse order that they were
    %plotted
    hln = hln(end:-1:1);
end

xl = get(hax,'XLim');
yl = get(hax,'YLim');

n = min(length(labels), length(hln));

x = [];
if (ischar(opts.location))
    switch lower(opts.location)
        case 'right'
            x = xl(2) * ones(size(hln));
            opts.HorizontalAlignment = 'left';
            
        case 'best'
            x = (xl(1) + 2/3 * (xl(2)-xl(1))) * ones(size(hln));
    end
elseif isnumeric(opts.location)
    if (numel(opts.location) == 1)
        x = opts.location * ones(size(hln));
    elseif (numel(opts.location) == numel(labels))
        x = opts.location;
    end
end

if (isempty(x))
    error('Unrecognized location parameter');
end

rot = [];
islinerot = false;
if (isnumeric(opts.rotation))
    if (numel(opts.rotation) == numel(hln))
        rot = opts.rotation;
    elseif (numel(opts.rotation) == 1)
        rot = opts.rotation * ones(size(hln));
    end
elseif ischar(opts.rotation)
    switch opts.rotation
        case 'line'
            islinerot = true;
            rot = zeros(size(hln));
    end
end
if (isempty(rot))
    error('Unrecognized rotation option');
end

AR = diff(yl) / diff(xl);

y = zeros(size(hln));
col = cell(size(hln));
htxt = -1 * ones(size(hln));
for i = 1:n
    xd = get(hln(i), 'XData');
    yd = get(hln(i), 'YData');
    
    [xd,ord] = sort(xd);
    yd = yd(ord);
    
    norep = [true diff(xd) > 0];
    xd = xd(norep);
    yd = yd(norep);
    
    if (x(i) > max(xd))
        x(i) = max(xd);
    elseif (x(i) < min(xd))
        x(i) = min(xd);
    end
    
    y(i) = interp1(xd,yd, x(i));
    col{i} = get(hln(i),'Color');
    
    htxt(i) = text(x(i),y(i), labels{i}, 'Color',col{i}, ...
        'BackgroundColor','w', 'EdgeColor',col{i}, ...
        'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
        'Parent',hax);
    
    if (islinerot)
        ext = get(htxt(i),'Extent');
        
        overind = last(xd < ext(1)):first(xd > ext(1)+ext(3));
        
        if ~isempty(overind)
            dy = mean(diff(yd(overind)));
            dx = range(xd(overind));
        
            rot(i) = atan2(dy, dx*AR)*180/pi;
        end
        
        set(htxt(i),'Rotation',rot(i));
    end
end

