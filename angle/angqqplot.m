function angqqplot(ang1,ang2)
% function angqqplot(ang)
%   or     angqqplot(ang1,ang2)
% Equivalent of a quantile-quantile plot for one or two angular data sets.
% For one data set, plots the data against the von Mises distribution.  
% For two datasets, plots the distributions (centered on their medians) 
% against one another.
%
% NB: this doesn't work right!
%
% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2010, Eric Tytell


ang1 = shiftdim(ang1(isfinite(ang1)));

if (nargin == 2),
  ang2 = shiftdim(ang2(isfinite(ang2)));
  if (length(ang1) > length(ang2)),
    swap = ang1;
    ang1 = ang2;
    ang2 = swap;
  end;

  z1 = sin(0.5*(ang1 - angmedian(ang1)));
  z2 = sin(0.5*(ang2 - angmedian(ang2)));

  z1 = sort(z1);
  z2 = sort(z2);

  k = round((1:length(ang1))*length(ang2)/length(ang1));
  m1 = min([z1; z2]);
  m2 = max([z1; z2]);
  plot([m1 m2],[m1 m2],'k:', z1,z2(k),'+');
elseif (nargin == 1),
  [mu] = angmean(ang1);
  kappa = angkappa(ang1);

  z1 = sin(0.5*(ang1 - mu));
  z1 = sort(z1);

  quant = (1:length(ang1))'/length(ang1);

  VMtheta = VMinv(quant, pi,kappa) - pi;
  zvm = sin(0.5*VMtheta);

  m1 = min([z1; zvm]);
  m2 = max([z1; zvm]);

  plot([m1 m2],[m1 m2],'k:', zvm,z1, '+');
end;


 