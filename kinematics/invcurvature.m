function [x,y] = invcurvature(curve,len,ang0)
% function [x,y] = invcurvature(curve,len,ang0)
% Estimates the midline coordinates based on curvature.  ang0 is optional, but it
% specifies the angle of the first body segment.  Otherwise the first angle is assumed to
% be zero.  len is the body length.

ds = len/(size(curve,1)-1);

dang = zeros(size(curve,1)-1,size(curve,2));
dang(2:end,:) = curve(2:end-1,:) * ds;
if (nargin == 3),
    dang(1,:) = ang0;
end;
segang = cumsum(dang);

dx = ds * cos(segang);
dy = ds * sin(segang);

x = zeros(size(curve));
y = zeros(size(curve));

x(2:end,:) = cumsum(dx);
y(2:end,:) = cumsum(dy);
