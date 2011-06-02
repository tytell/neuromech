function [burstctr,burstind] = findbursts_timing(t,spike,opt)

%check whether they specified things in terms of frequency
if (isfield(opt,'burstfreq')),
    %baseline burst duration and interburst duration            
    %minimum cutoffs will be a fraction of those
    opt.cycledur = 1/burstfreq;
    if (~isfield(opt,'burstdur')),
        burstdur = opt.cycledur * opt.dutycyle;
        burstdur = burstdur * opt.cutoff;
    end;
    if (~isfield(opt,'interburstdur')),
        interburstdur = opt.cycledur - burstdur;
        interburstdur = interburstdur * opt.cutoff;
    end;
end;
if (isfield(opt,'cycledur') || isempty(opt.cycledur)),
    opt.cycledur = burstdur + interburstdur;
end;

if (strcmp(opt.method,'timing')),
    sig = diff(t);
    thresh = interburstdur;
else
    sig = 1./spike;
    if (isempty(opt.threshold)),
        opt.threshold = getthreshold(t,sig);
    end;
    opt.minspikes = 0;
    thresh = 1./opt.threshold;
end;

%find regions where something (spike rate, or integrated signal) is
%more than a threshold
up = find((sig(1:end-1) > thresh) & ...
          (sig(2:end) <= thresh)) + 1;
down = find((sig(1:end-1) <= thresh) & ...
            (sig(2:end) > thresh)) + 1;

if (isempty(up) || isempty(down)),
    warning('No bursts found.');
    varargout = {[],zeros(2,0)};
    return;
end;

%discard extra ups or downs at the beginning or end of the sequence
if (down(1) < up(1)),
    down = down(2:end);
end;
if (up(end) > down(end)),
    up = up(1:end-1);
end;

%merge bursts that are separated by less than mergegap
merge = t(up(2:end))-t(down(1:end-1)) < opt.mergegap;
up = up([true ~merge]);
down = down([~merge true]);

%check to make sure the bursts are long enough and contain enough
%spikes
longenough = ((down - up) >= opt.minspikes) & ...
    (t(down) - t(up)) >= burstdur;
up = up(longenough);
down = down(longenough);

%get the burst center positions
ctr = zeros(size(up));
for i = 1:length(up),
    ctr(i) = mean(t(up(i):down(i)));
end;

if (opt.correctmissedbursts),
    medper = nanmedian(diff(ctr));
end;
addup = [];
adddown = [];
addctr = [];

%now eliminate multiple bursts within the minimum cycle period
nbursts = length(up);

good = true(size(up));
i = 1;
while (i <= length(up)-1),
    nextfew = i+1:i+10;
    nextfew = nextfew(nextfew <= nbursts);
    
    withinper = find(good(nextfew) & (ctr(nextfew) - ctr(i) < opt.cycledur));
    if (isempty(withinper)),
        if (opt.correctmissedbursts),
            per = ctr(i+1) - ctr(i);
            nskip = round(per / medper);
            if (nskip == 2),
                k = down(i)+1:up(i+1)-1;

                up1 = find((sig(k-1) > thresh/2) & ...
                           (sig(k) <= thresh/2)) + down(i);
                down1 = find((sig(k) <= thresh/2) & ...
                             (sig(k+1) > thresh/2)) + down(i);
                
                if (~isempty(up1) && (length(up1) == length(down1))),
                    ctr1 = zeros(size(up1));
                    for j = 1:length(up1),
                        ctr1(j) = mean(t(up1(j):down1(j)));
                    end;

                    %choose the one closest to where the burst
                    %ought to be
                    expected = (ctr(i+1) + ctr(i))/2;
                    [q,best] = min(abs(ctr1 - expected));
                    
                    addup(end+1) = up1(best);
                    adddown(end+1) = down1(best);
                    addctr(end+1) = ctr1(best);
                end;
            end;
        end;
        i = i+1;
    else
        nspike = down(i+withinper) - up(i+withinper);
        [m,ind] = max(nspike);
        good(i+withinper) = false;
        good(i+ind) = true;

        i = i+ind;
    end;
end;

up = sort([up(good) addup]);
down = sort([down(good) adddown]);
burstctr = sort([ctr(good) addctr]);

burstind = [up; down];


