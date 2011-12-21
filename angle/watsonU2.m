function P = watsonU2(U2,n1,n2)

U = load('WatsonU2.mat');

n = sort([n1 n2]);

if (n(1) < 4),
    error('Cannot calculate critical U2 for n < 4.');
end;

loga = log(U.alpha);

if (n(1) <= 12),
    i = find(U.n1 == n(1));
    j = find(U.n2(i) <= n(2));
    j = j(end);
elseif (n(1) > 100),
    i = size(U.U2crit,1);
    j = 1;
else
    [m,i] = min(abs(n(1) - U.n1(1)));
    j = 1;
end;

if (U2 < U.U2crit(i(j),1)),
    P = 1;
elseif (U2 > max(U.U2crit(i(j),:))),
    P = 0;
else
    logP = interp1(U.U2crit(i(j),:),loga,U2);
    P = exp(logP);
end;
