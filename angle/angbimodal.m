function [mu1,mu2, k1,k2, p] = angbimodal(ang)

ang = ang(isfinite(ang));
ang = shiftdim(ang);

mu1 = angmean(ang);
mu2 = mod(mu1 + pi, 2*pi);
k1 = angkappa(ang);
k2 = k1;
p = 0.5;

n = length(ang);
for i = 1:3,
  C(i) = sum(cos(i*ang)/n);
  S(i) = sum(sin(i*ang)/n);
end;

x0 = [mu1 mu2 k1 k2 p];
x = fminsearch(@VMbimodal,x0,[],C,S);

mu1 = x(1);
mu2 = x(2);
k1 = x(3);
k2 = x(4);
p = x(5);

if (k1 < 0),
  mu1 = mod(mu1 + pi, 2*pi);
  k1 = -k1;
end;
if (k2 < 0),
  mu2 = mod(mu2 + pi, 2*pi);
  k2 = -k2;
end;

function r = VMbimodal(x, C,S)

mu1 = x(1);
mu2 = x(2);
k1 = x(3);
k2 = x(4);
p = x(5);

for i = 1:3,
  dC(i) = C(i) - (p*A(i,k1)*cos(i*mu1) + (1-p)*A(i,k2)*cos(i*mu2));
  dS(i) = S(i) - (p*A(i,k1)*sin(i*mu1) + (1-p)*A(i,k2)*sin(i*mu2));
end;

r = sum(dC.^2) + sum(dS.^2);

function a = A(i,k)

I0 = besseli(0,k);
I1 = besseli(1,k);

A1 = I1/I0;
a = A1;
if (i > 1),
  A2 = 1 - 2/k*A1;
  a = A2;
  if (i > 2),
    A3 = A1 - 4/k*A2;
    a = A3;
  end;
end;


