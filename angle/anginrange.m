function inrange = anginrange(a, lo,hi)

if (mod(lo - hi,2*pi) == pi)
    warning('Results may be undefined for range = pi');
end

Clo = cos(lo);
Slo = sin(lo);

Chi = cos(hi);
Shi = sin(hi);

C = cos(a);
S = sin(a);

locrosshi = Clo .* Shi - Slo .* Chi;
locrossa = Clo .* S - Slo .* C;
acrosshi = C .* Shi - S .* Chi;

inrange = (sign(locrosshi) == sign(locrossa)) & (sign(locrosshi) == sign(acrosshi));
