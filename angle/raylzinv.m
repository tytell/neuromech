function P = raylzinv(Z,n)

if (any(n(:) < 10)),
  warning('Accuracy will be low with low n.');
end;

R = sqrt(Z.*n);
P = exp(sqrt(1+4*n+4*(n.^2-R.^2)) - (1+2*n)); % from Zar, p 617


