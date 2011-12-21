function s = angdev(ang,dim)
% function s = angdev(ang,dim)
% Calculates sample angular standard deviation
% From Fisher 1993, page 32

if (nargin == 1),
  dim = 1;
end;

angx = nanmean(cos(ang),dim);
angy = nanmean(sin(ang),dim);

r = sqrt(angx.^2 + angy.^2);
lr = log(r);
lr(lr > 0) = NaN;

s = sqrt(-2*lr);
