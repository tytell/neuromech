function [P,F,df1,df2] = angpairedtest(ang1,ang2)
% Performs a paired test for diffences in two vectors of angles (in
% radians).  Each element in ang1 must have a corresponding element in
% ang2.
% Algorithm from Zar 1999, p. 645

angx = cos(ang1) - cos(ang2);
angy = sin(ang1) - sin(ang2);

nonan = find(isfinite(angx) & isfinite(angy));
angx = angx(nonan);
angy = angy(nonan);

k = sum(isfinite(angx));

sx2 = sum(angx.^2) - sum(angx).^2/k;
sy2 = sum(angy.^2) - sum(angy).^2/k;
sxy = sum(angx.*angy) - sum(angx).*sum(angy)/k;

xbar = mean(angx);
ybar = mean(angy);

F = k*(k-2)/2*(xbar^2*sy2 - 2*xbar*ybar*sxy + ybar^2*sx2) / ...
    (sx2*sy2 - sxy^2);

df1 = 2;
df2 = k-2;
P = 1 - fcdf(F,df1,df2);
