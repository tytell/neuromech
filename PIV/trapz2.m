function f = trapz2(x,y,z)
% function f = trapz2(x,y,z)
% 2D trapezoidal integration, dealing with NaNs as if they bounded the
% integration.  Works best when z defines a convex region bounded by NaNs.

good = isfinite(z);
r = size(good,1);
c = size(good,2);
i = 2:r-1;
j = 2:c-1;

%a block with linear edges defined at the corners by f11, f12, f21, f22
%will have a curved surface and a volume of 
%0.25*h^2*(f11 + f12 + f21 + f22), where h is the xy distance between the
%function values.  If we sum this over a bunch of adjacent points, we get
%the 2D integral.  Interior points are summed more than once.

%figure out how many non-NaN neighbors a point has.  This will define how
%many times the point would be repeated
nbh = zeros(size(z));
nbh(i,j) = good(i-1,j) + good(i+1,j) + good(i,j-1) + good(i,j+1);
nbh([1 r],j) = good([2 r-1],j) + good([1 r],j-1) + good([1 r],j+1);
nbh(i,[1 c]) = good(i-1,[1 c]) + good(i+1,[1 c]) + good(i,[2 c-1]);
nbh([1 r],[1 c]) = good([2 r-1],[1 c]) + good([1 r],[2 c-1]);

%weights
w = zeros(size(z));
w(nbh == 2) = 0.25;
w(nbh == 3) = 0.5;
w(nbh == 4) = 1;

%NB: we ignore solitary points (nbh == 0) and isolated groups of two (nbh
%== 1).  Groups with 3 give the corner a weight of 0.25 and ignore the
%others.  Probably this is incorrect, but it would be difficult and
%computationally intensive to identify the bad points.

%do the sum
dx = x(1,2)-x(1);
dy = y(2,1)-y(1);
f = dx*dy * sum(w(good) .* z(good));
