function [indpeak, per,amp] = analyzeTailbeat(fishlen,t, hx,hy, tx,ty, varargin)

opt.smoothper = 1;      % in sec
opt.smoothmethod = 'loess';
opt.ampthresh = 0.01;   % fraction of body length
opt.npt = 20;

opt = parsevarargin(opt,varargin, 7);

dt = t(2)-t(1);
assert(~isnan(dt));

ang = atan2(ty-hy, tx-hx);
if (range(ang) > pi)
    %check if the angle wraps around
    ang = unwrap(ang);
end

%smooth the angle to get the average body angle
good = isfinite(ang);
bodyang = NaN(size(hx));

smoothfrac = opt.smoothper / (sum(good)*dt);
if (smoothfrac < 1)
    bodyang(good) = smooth(ang(good),smoothfrac,opt.smoothmethod);
else
    p = polyfit(t(good),ang(good),2);
    bodyang(good) = polyval(p,t(good));
end

tailang = ang - bodyang;

thresh = atan(opt.ampthresh);

[phase,cycleind,cyclet] = getphase(t,tailang,'threshold',thresh);

good = cycleind ~= 0;
if all(~good)
    warning('No tail beats detected.');
    indpeak = NaN(opt.npt,1);
    amp = NaN(opt.npt,1);
    per = NaN(opt.npt,1);
else
    angamp = NaN(size(cycleind));
    angamp(good) = tailang(cycleind(good));

    goodcycle = any(abs(angamp) > thresh);
    cycleind = cycleind(:,goodcycle);
    cyclet = cyclet(:,goodcycle);

    per1 = NaN(size(cyclet));
    per1(2,:) = cyclet(3,:) - cyclet(1,:);
    per1(4,1:end-1) = cyclet(1,2:end) - cyclet(3,1:end-1);

    good = cycleind ~= 0;
    txcycle = NaN(size(cycleind));
    txcycle(good) = tx(cycleind(good));
    tycycle = NaN(size(cycleind));
    tycycle(good) = ty(cycleind(good));

    amp1 = NaN(size(cyclet));
    amp1(2,:) = 0.5*(sqrt((txcycle(2,:)-txcycle(1,:)).^2 + (tycycle(2,:)-tycycle(1,:)).^2) + ...
        sqrt((txcycle(2,:)-txcycle(3,:)).^2 + (tycycle(2,:)-tycycle(3,:)).^2));
    amp1(4,1:end-1) = 0.5*(sqrt((txcycle(4,1:end-1)-txcycle(3,1:end-1)).^2 + ...
        (tycycle(4,1:end-1)-tycycle(3,1:end-1)).^2) + ...
        sqrt((txcycle(4,1:end-1)-txcycle(1,2:end)).^2 + ...
        (tycycle(4,1:end-1)-tycycle(1,2:end)).^2));

    indpeak = NaN(opt.npt,2*size(cycleind,2));
    indpeak(end,:) = flatten(cycleind([2 4],:));
    indpeak(indpeak == 0) = NaN;

    per = NaN(size(indpeak));
    per(end,:) = flatten(per1([2 4],:));

    amp = NaN(size(indpeak));
    amp(end,:) = flatten(amp1([2 4],:));

    good = any(isfinite(indpeak));
    indpeak = indpeak(:,good);
    amp = amp(:,good);
    per = per(:,good);
end


