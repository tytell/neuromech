function [I1,I2,xp,yp,up,vp] = imparticle(x,y,u,v,noise, w,h, ...
                                          ppp,d,stdd, pint,pfreq,stdint, mode)
% function [I1,I2,xp,yp,up,vp] = imparticle(x,y,u,v,noise, w,h, ...
%                                     ppp,d,stdd, pint,pfreq,stdint, mode)
% Builds two artificial PIV images based on the velocity field (x,y,u,v).
% Generates randomly positioned particles with gaussian intensity
% distributions and moves them according to the specified velocity field.
%
% Adds normally distributed noise to the velocity vectors with standard
% deviation noise.  The returned images have size [h w].  The particles are
% distributed with density ppp (in particles per pixel), have diameter d,
% and standard deviation of diameters stdd.  They also have a range of
% intensities pint, each with frequency pfreq.  Both pint and pfreq should
% be vectors of the same length.  Particle intensities in the first image
% vary with standard deviation stdint.  mode determines how the particle
% intensities are calculated.  mode = 1 simply takes the value of the
% Gaussian intensity at the center of the pixel for each pixel (fastest).
% mode = 2 integrates the Gaussian intensity over each pixel (most correct,
% but slower).  Mode 2 is particularly important when particle diameters are
% less than 2 pixels.
%
% Additionally returns the positions and velocities of all the particles.
% This can be useful if you want to generate multiple images of the same set
% of particles moving.

if (nargin == 13)
  mode = 1;
end;

if (any(isnan(x(:))) | any(isnan(y(:))) | any(isnan(u(:)) | isnan(v(:)))),
  error('Cannot construct spline with NaNs.');
end;

tic;

d = shiftdim(d);
if (prod(size(d)) == 1),
  d = [d; d];
end;
stdd = shiftdim(stdd);
if (prod(size(stdd)) == 1),
  stdd = [stdd; stdd];
end;

partsize = round(d./2+stdd./2);    % size of an individual particle image
partsize(partsize == 0) = 1;

if ((ndims(u) == 2) & (any(size(u) == 1))),		% vector velocity
  k = find((x >= partsize(1)+1) & (x <= w-partsize(1)) & ...
           (y >= partsize(2)+1) & (y <= h-partsize(2)));
  np = length(k);
  
  xp = x(k);
  yp = y(k);
  up = u(k);
  vp = v(k);

  if (prod(size(d)) == 2),
    dp = repmat(d,[1 np]);
  else
    dp = d(:,k);
  end;
else,
  x1 = max(min(x(:)), partsize(1)+1);
  y1 = max(min(y(:)), partsize(2)+1);
  x2 = min(max(x(:)), w-partsize(1));
  y2 = min(max(y(:)), h-partsize(2));
  
  if (~any(size(x) == 1)),		% x is not a vector
    x = x(1,:);
    y = y(:,1)';
  end;
  
% make randomly arranged particles
  np = round(w*h*ppp);
  xp = rand(1,np)*(x2-x1) + x1;
  yp = rand(1,np)*(y2-y1) + y1;

  datap = prod(size(u))/(range(x)*range(y));
  if (datap < ppp/5),
    % data density is low -- use a good spline function to produce nice
    % vectors
    
    % note transposition of u and v -> spline functions use (x,y) referencing
    % not (i,j)
    U = zeros(2,size(u,2),size(u,1));
    U(1,:,:) = u';
    U(2,:,:) = v';
    U(isnan(U)) = 0;
    
    sp = csapi({x, y}, U);
    Ur = fnval(sp, [xp; yp]);
    
    up = squeeze(Ur(1,:)) + randn(1,np)*noise;
    vp = squeeze(Ur(2,:)) + randn(1,np)*noise;
  else,
    % data density is high -- use a table interpolation routine so that we can
    % be fast
    
    up = interp2(x,y,u, xp,yp, '*cubic') + randn(1,np)*noise;
    vp = interp2(x,y,v, xp,yp, '*cubic') + randn(1,np)*noise;
  end;
  
  % calculate particle sizes
  dp = repmat(d,[1 np]) + repmat(randn(1,np),[2 1]).*repmat(stdd,[1 np]);
end;

% sort out the particle intensities
pfreq = ceil(pfreq*np/sum(pfreq));
p = 0;
for i = 1:length(pfreq),
  partint(1,p+(1:pfreq(i))) = pint(i);
  p = p+pfreq(i);
end;
%take the first np intensities (shortchanging the last one by a few
%particles sometimes), because the rounding doesn't always add up right
partint = partint(1:np);

I1 = makepartim(w,h, xp,yp, partint, partsize,dp, mode);

xp2 = xp + up;
yp2 = yp + vp;
partint = partint + randn(size(partint))*stdint;
partint(partint > 1) = 1;
k = find((xp2 >= partsize(1)+1) & (xp2 <= w-partsize(1)) & ...
         (yp2 >= partsize(2)+1) & (yp2 <= h-partsize(2)) & ...
         (partint > 0));
% fprintf(1, '%f%% (%d) particles lost\n', (np - length(k))/np*100, (np-length(k)));

I2 = makepartim(w,h, xp2(k),yp2(k), partint(k), partsize,dp(:,k), mode);

deltat = toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% makepartim - Makes the image of particles at (xp,yp) with intensities
% partint, pixel size partsize, gaussian size d
% This is a hard problem.  Each particle is a small image (partsize x
% partsize) that needs to be added to the overall image at a random place.
% Since the number of particles np tends to be large, referencing that many
% random locations is very time consuming.  So we do a more complicated
% method: Each particle image consists of npixpart pixels that have certain
% indices around the center of the particle itself (xp,yp).  Thus, we know
% the indices within the overall image that need to be set to the various
% intensities of the particle images.  So we generate the particle images in
% a linear fashion to be a (npixpart x np) matrix, with locations in the
% image xpart and ypart Then we generate the columnwise indices
% corresponding to all of xpart and ypart.  Ideally, we could just then set
% the image at those indices to partim.  But, particles often overlap,
% meaning that the indices repeat.  So we have to deal with repeated
% indices, which (mostly) means summing up the particle intensities for each
% repeated index and assigning the image once at that index to the sum of
% the appropriate intensities.  It's complicated, but see below.
function I = makepartim(w,h,xp,yp,partint, partsize,d, mode)

I = zeros(h,w);

[xpart,ypart] = meshgrid(-partsize(1):(partsize(1)+1), ...
                         -partsize(2):(partsize(2)+1));
npixpart = prod(size(xpart));        % number of pixels

% integer pixel coordinates
ixp = floor(xp);
iyp = floor(yp);
% fractional displacements
dxp = xp-ixp;
dyp = yp-iyp;

% for each particle, generate a column with all the indices of the pixels
% it will be defined on, relative to its center
xpart = repmat(xpart(:),[1 length(xp)]);
ypart = repmat(ypart(:),[1 length(yp)]);

% then define the image at those pixels, with the fractional displacement
partint0 = 1/4*(erf(sqrt(8).*(-0.5)./d(1)) - erf(sqrt(8).*(+0.5)./d(1))) .* ...
    (erf(sqrt(8).*(-0.5)./d(1)) - erf(sqrt(8).*(+0.5)./d(1)));

partintl = repmat(partint,[npixpart 1]);
xx = xpart - repmat(dxp,[npixpart 1]);
yy = ypart - repmat(dyp,[npixpart 1]);
ddx = repmat(d(1,:),[npixpart 1]);
ddy = repmat(d(2,:),[npixpart 1]);
if (mode == 1),
  partim = exp(-8*(xx.^2./ddx.^2 + yy.^2./ddy.^2));
elseif (mode == 2),
  partim = 1/(4*partint0) * ...
           (erf(sqrt(8).*(xx-0.5)./ddx) - erf(sqrt(8).*(xx+0.5)./ddx)) .* ...
           (erf(sqrt(8).*(yy-0.5)./ddy) - erf(sqrt(8).*(yy+0.5)./ddy));
end;
partim = partim .* partintl;

% add on the integer coordinates
xpart = xpart + repmat(ixp,[npixpart 1]);
ypart = ypart + repmat(iyp,[npixpart 1]);

% remove values effectively equal to zero
q = find(partim > 1/256);

k = ypart(q) + (xpart(q)-1)*h;
partim = partim(q);
I(k) = partim;

% some particles may overlap, so we have to deal with them
% get the indices of the pixel position for each particle and see which are
% the same
[k,isort] = sort(k);
partim = partim(isort);

good = [1; find(diff(k) ~= 0)+1];     % non repeated
I(k(good)) = partim(good);

% find all the repeats
rep = find(diff(k) == 0)+1;

if (~isempty(rep)),
  % sum up the repeated indices
  repsum = cumsum(partim(rep));
  % find the points where repeats stop
  % i.e. k = [1 1 2 2 3 4 5 5 6] would eventually give endrep = [2 4 5 6 8]
  endrep = find(diff(rep) > 1);
  % the difference of the cumsum at each of the endrep points gives the
  % sum of all the repeated points for each pixel
  repsum = diff([0; repsum([endrep; end])]);
  % set the image appropriately
  I(k(rep([endrep; end]))) = repsum + I(k(rep([endrep; end])));
end;

% elseif (mode == 2),
% % step through single repeats, then double and so forth, each time adding
% % on the next value
% 	while (~isempty(rep)),
%         repval = k(rep);
%         
%         good = [1; find(diff(repval) ~= 0)+1];
%         I(repval(good)) = I(repval(good)) + partim(rep(good));
%         
%         rep = rep(diff(repval) == 0) + 1;
% 	end;
% end;

I(I > 1) = 1;

