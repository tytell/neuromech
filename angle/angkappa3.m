function k = angkappa3(u,v,w)

p = 3;                                  % number of dimensions

sz = size(u);
k = zeros([1 sz(2:end)]);

R = sqrt(mean(u).^2 + mean(v).^2 + mean(w).^2);

big = find(R > 0.9);
k(big) = (p-1)./(2*(1-R(big)));

small = find(R < 0.05);
k(small) = p*R(small).*(1 + p/(p+2)*R(small).^2 + ...
                        p^2*(p+8)/(p+2)^2*(p+4)*R(small).^4);

% spline based inverse

ksmall = p*0.05.*(1 + p/(p+2)*0.05.^2 + ...
                        p^2*(p+8)/(p+2)^2*(p+4)*0.05.^4);
kbig = (p-1)./(2*(1-0.9));

k0 = ksmall:0.1:kbig;
if (p == 3),
    Ap = coth(k0) - 1./k0;

    rest = find((R >= 0.05) & (R <= 0.9));
    k(rest) = spline(Ap,k0, R(rest));
end;

