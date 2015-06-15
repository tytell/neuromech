function isdivisible = divisible(x,y, tol)

if nargin == 2
    tol = 1e-10;
end

q = x/y;
isdivisible = abs(q - round(q)) < tol;
