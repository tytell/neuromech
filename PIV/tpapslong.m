function [st,p] = tpapslong(x,y,p)
%TPAPS Thin-plate smoothing spline.
%
%   F = TPAPS(X,Y)  is the stform of a thin-plate smoothing spline  f  for
%   the given data sites X(:,i) and corresponding data values Y(:,i). 
%   The X(:,i) must be distinct points in the plane, and X and Y must have
%   the same number of columns. 
%   The thin-plate smoothing spline  f  is the unique minimizer of the 
%   weighted sum
%
%                     P*E(f) + (1-P)*R(f) ,
%
%   with E(f) the error measure
%
%       E(f) :=  sum_j { | Y(:,j) - f(X(:,j)) |^2 : j=1,...,n }
%
%   and R(f) the roughness measure
%
%       R(f) := integral  (D_1 D_1 f)^2 + 2(D_1 D_2 f)^2 + (D_2 D_2 f)^2.
%
%   Here, the integral is taken over the entire 2-space, and 
%   D_i denotes differentiation with respect to the i-th argument, hence 
%   the integral involves the second derivatives of  f .
%   The smoothing parameter P is chosen in an ad hoc fashion
%   in dependence on the sites X.
%
%   TPAPS(X,Y,P) provides the smoothing parameter P, a number expected to 
%   be between  0  and  1 .  As P varies from 0 to 1, the smoothing spline
%   changes, from the least-squares approximation to the data by a linear
%   polynomial when P is 0, to the thin-plate spline interpolant to the data
%   when P is 1.
%
%   [F,P] = TPAPS(...) also returns the smoothing parameter used.
%
%   Warning: The determination of the smoothing spline involves the solution 
%   of a linear system with as many unknowns as there are data points.
%   Since the matrix of this linear system is full, the solving can take a long
%   time even if, as is the case here, an iterative scheme is used when there
%   are more than 728 data points. The convergence speed of that iteration is
%   strongly influenced by P, and is slower the larger P is. So, for large
%   problems, use interpolation (P equal to 1) only if you can afford the time.
%
%   Examples:
%
%      rand('seed',23); nxy = 31;
%      xy = 2*(rand(2,nxy)-.5); vals = sum(xy.^2);
%      noisyvals = vals + (rand(size(vals))-.5)/5;
%      st = tpaps(xy,noisyvals); fnplt(st), hold on
%      avals = fnval(st,xy);
%      plot3(xy(1,:),xy(2,:),vals,'wo','markerfacecolor','k')
%      quiver3(xy(1,:),xy(2,:),avals,zeros(1,nxy),zeros(1,nxy), ...
%               noisyvals-avals,'r'), hold off
%   generates the value of a very smooth function at 31 random sites,
%   adds some noise to it, then constructs the smoothing spline to these
%   noisy data, plots the smoothing spline, the exact values (as black
%   balls) the smoothing is trying to recover, and the arrow leading from
%   the smoothed values to the noisy values.
%
%      n = 64; t = linspace(0,2*pi,n+1); t(end) = [];
%      values = [cos(t); sin(t)];
%      centers = values./repmat(max(abs(values)),2,1);
%      st = tpaps(centers, values, 1);
%      fnplt(st), axis equal
%   constructs a map from the plane to the plane that carries the unit square,
%   {x in R^2: |x(j)|<=1, j=1:2}, pretty much onto the unit disk 
%   {x in R^2: norm(x)<=1}, as shown by the picture generated.
%
%   See also CSAPS, SPAPS.

%   Copyright 1987-2001 C. de Boor and The MathWorks, Inc.
%   $Revision$  $Date$

[d,nx] = size(x);
dd = d+1;
if nx<dd
   error(['Thin-plate spline smoothing requires at least ',num2str(dd), ...
          ' data sites.'])
end

[dy,ny] = size(y);
if ny~=nx, error(' x  and  y  must have the same number of columns.'), end

[Q,R] = qr([ x.' ones(nx,1)]);
radiags = sort(abs(diag(R)));
if radiags(1)<1.e-14*radiags(end)
   error('Some nontrivial linear polynomial vanishes at all data sites.')
end

if nx==3 % simply return the interpolating plane
   st = stmak(x,[zeros(dy,3), y/(R(1:dd,1:dd).')],'tp00');
   p = 1;
   return
end

Q1 = Q(:,1:dd); Q(:,1:dd) = [];
if nargin==3&~isempty(p)&p==0 % get the linear least squares polynomial:
   st = stmak(x, [zeros(dy,nx), (y*Q1)/(R(1:dd,1:dd).')],'tp00');
   return
end

colmat = stcol(x,x,'tr');

if nargin<3|isempty(p) % we must supply the smoothing parameter
  p = 1/(1+mean(diag(Q'*colmat*Q)));
end

if p~=1, colmat(1:nx+1:nx^2) = colmat(1:nx+1:nx^2)+(1-p)/p; end
coefs1 = (y*Q/(Q'*colmat*Q))*Q';
coefs2 = ((y - coefs1*colmat)*Q1)/(R(1:dd,1:dd).');

st = stmak(x,[coefs1,coefs2],'tp00');


function vals = tppval(x,st,Q2,p)
%TPPVAL evaluation for iterative solution of thin-plate spline smoothing system

if isempty(Q2)
   st.coefs = x.';
   vals = stval(st,st.centers).';
else
   st.coefs = (Q2*x).';
   vals = (stval(st,st.centers)*Q2).';
end
if p~=1 % TPPVAL is never called when p==0
   vals = vals + x*((1-p)/p);
end

