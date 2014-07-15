function dphase = angdiff(phase1,phase2)

dphase = phase1 - phase2;
if (dphase < 0)
    dphase2 = 2*pi + dphase;
else
    dphase2 = dphase - 2*pi;
end

isflip = abs(dphase2) < abs(dphase);
dphase(isflip) = dphase2(isflip);
