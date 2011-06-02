function widthL = digitizeFishWidth(I, n)
% function widthL = digitizeFishWidth(I, n)
% Digitizes the fish width from an image I, producing n evenly spaced
% points.
%
% Mercurial revision hash: $Revision: dd099ded8eec $ $Date: 2010/08/11 13:58:29 $
% Copyright (c) 2010, Eric Tytell

if (nargin == 1),
    n = 20;
end;

imshow6(I,'n');
fprintf('Zoom to fish body and press return.');
zoom on;
pause;

zoom off;
done = false;
while ~done,
    fprintf('Click left side of fish body\n');
    [lx,ly] = ginputb;

    fprintf('Click right side of fish body\n');
    [rx,ry] = ginputb;

    %calculate arc length on each side
    ls = [0 cumsum(sqrt(diff(lx).^2 + diff(ly).^2))];
    rs = [0 cumsum(sqrt(diff(rx).^2 + diff(ry).^2))];

    %interpolate same points along the arc length
    rx2 = spline(rs,rx, ls);
    ry2 = spline(rs,ry, ls);

    %now estimate the midline
    mx0 = (lx + rx2)/2;
    my0 = (ly + ry2)/2;

    %spline will give NaNs if it tries to extrapolate, so make sure the
    %final point is OK
    if (~isfinite(mx0(end)) || ~isfinite(my0(end))),
        mx0(end) = (lx(end) + rx(end))/2;
        my0(end) = (ly(end) + ry(end))/2;
    end;

    %and interpolate so that they're equally spaced along the arc
    ms0 = [0 cumsum(sqrt(diff(mx0).^2 + diff(my0).^2))];
    ms = linspace(0,ms0(end),n);
    fishlen = ms(end);

    width0 = sqrt((lx - rx2).^2 + (ly - ry2).^2);
    widthL = spline(ms0, width0/fishlen, ms);

    %now plot it -- assuming the fish is straight
    ms1 = linspace(0,1,n);
    plot(ms1,widthL/2, ms1,-widthL/2);
    axis equal;

    done = inputyn('Is outline OK? ');
end;
