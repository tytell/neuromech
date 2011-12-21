function P = ang3rayltest(R,N)

N = sum(isfinite(angx) & isfinite(angy) & isfinite(angz));
R = sqrt(nansum(angx).^2 + nansum(angy).^2 + nansum(angz).^2);

if (N < 4),
    warning('Probabilities will be erroneous for low N');
elseif (N >= 100),
    P = 1 - chi2cdf(3*R^2/N,2);
else
    P = quad(@rayldist,R,N,[],[],N);
end;



function f = rayldist(R, N)

n = length(R);
s = repmat((0:floor(N/2))',[1 n]);
R = repmat(R,[size(s,1) 1]);

xbrak = N - R - 2*s;
xbrak(xbrak < 0) = 0;
ff = combination(N,s) .* (-1).^s .* xbrak.^(N-2);

f = R(1,:)./(2^(N-1)*factorial(N-2)) .* sum(ff);


