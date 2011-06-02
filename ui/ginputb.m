function [x,y,h] = ginputb(n, style)
% function [x,y,h] = ginputb(n, style)
% Like ginput, but beeps and plots a point at each click
% n is the number of points to collect, style is the point style
% to display (default is 'r+')
% Both style and n are optional.
%
% [x,y] = ginput
%   Collects points until return is pressed
%
% Returns the points in x,y and handles to the graphics objects
% in h.  That way you can get rid of them without replotting the
% image (use delete(h))
%
% Mercurial revision hash: $Revision: dd099ded8eec $ $Date: 2010/08/11 13:58:29 $
% Copyright (c) 2010, Eric Tytell

% Sort out the number of parameters
if (nargin < 2)
    style = 'r+';
    if (nargin == 0)
        n = Inf;
    end;
end;

% If we're looking for a definite number of points, set
% up the arrays in advance.  It's quicker that way.
if (isfinite(n))
    x = zeros(n,1);
    y = zeros(n,1);
end;

% Collect the points
i = 1;

while (i <= n),
    [x1,y1] = ginput(1);
    if (isempty(x1))
        break;          % break out of the loop if they press return
    end;
    beep;
    x(i) = x1;
    y(i) = y1;
    if (i == 1),
        h = line('XData',x,'YData',y, 'MarkerEdgeColor',style(1), ...
            'Marker',style(2));
    else
        set(h,'XData',x,'YData',y);
    end;
    i = i+1;
end;

% Only return the number of points they actually clicked.  May
% be less than n.
x = x(1:i-1);
y = y(1:i-1);

