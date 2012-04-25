function P = ang3rayltest(R,N)
% ANG3RAYLTEST - Rayleigh test for significance of 3D angles
%  Doesn't work

% Mercurial revision hash: $Revision$ $Date$
% See http://rcn.ccs.tulane.edu/index.php5/Tytell_Matlab
% Copyright (c) 2012, Eric Tytell <tytell at jhu dot edu>

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


