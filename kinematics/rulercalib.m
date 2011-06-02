function [scale,err,N,distpix,distmm] = rulercalib(im,space)
% function [scale,err,N,distpix,distmm] = rulercalib(im,space)
% Calibrates an image based on points clicked on a ruler.  Allows the user
% to zoom to the ruler and then click as many points as he/she wants.  Uses
% the distances between the points to get the scale.
% Takes the image im and the distance (in mm) between the points clicked on
% the ruler.
% Returns the scale (in mm/pix), the error on the scale, the number of points
% used to calculate the scale, and the actual distances between points on the
% ruler and the expected distances in mm.
% NB: assumes that pixels are square.  The ruler can be at any angle in the frame,
% but it doesn't allow any pixel aspect ratio.
%
% Mercurial revision hash: $Revision: 18f43cd9074e $ $Date: 2010/08/10 21:11:58 $
% Copyright (c) 2010, Eric Tytell

if (ischar(im)),
	im = imread(im);
end;

imshow6(im,'n');
input('Zoom to ruler and press return');

fprintf(1,'Click on %d mm spaced marks',space);
[sx,sy] = ginput;
n = length(sx);

distmm = (0:n-1)'*space;
distpix = sqrt((sx - sx(1)).^2 + (sy - sy(1)).^2);

%fit the data (avoiding the statistics toolbox) and get the error.  We assume that the
%fit passes through zero, so we don't include a column of ones in the X parameter
[b,b_err] = lscov(distpix,distmm);
scale = b(1);
err = b_err(1);
N = n;
