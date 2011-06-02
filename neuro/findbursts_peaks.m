function [burstctr,burstind,peakind,peakrate,tbin,spikerate,ismerged] = findbursts_peaks(t,spike, opt)

%get the spike rate
if (numel(opt.binsize) > 1),
    binsize = opt.binsize(2)-opt.binsize(1);
else
    binsize = opt.binsize;
end;

[tbin,spikerate] = firingrate(t,opt.binsize,0);
smoothrate = runavg(spikerate,opt.smooth,2);

%sort out options to findpeaks
findpeakopt = {};
if (isfield(opt,'threshold') && ~isempty(opt.threshold)),
    findpeakopt = [findpeakopt,{'thresh',opt.threshold}];
else
    opt.threshold = 0;
end;
if (isfield(opt,'maxflat') && ~isempty(opt.maxflat)),
    findpeakopt = [findpeakopt,{'maxflat',opt.maxflat}];
end;
if (isfield(opt,'interburstdur') && ~isempty(opt.interburstdur)),
    %if peaks are closer than the interburstdur, then we *really* don't have separate
    %bursts.  But really, we need to check if the offset of one burst is closer to the
    %onset of the next than interburstdur, and findpeaks2 can't do that
    findpeakopt = [findpeakopt,{'minpeakdistance',...
                   opt.interburstdur / (tbin(2)-tbin(1))}];
end;

[peakrate,peakind] = findpeaks2(smoothrate,'max',findpeakopt{:});
peakind = peakind(:)';
peakrate = peakrate(:)';

%this gives us the index of the first spike in each bin
spikeind = round(cumsum(spikerate*binsize) + 1 - spikerate*binsize/2);

%and this gives us the index of the first spike in the bin that corresponds to the peak
%firing rate in each burst
peakspikeind = spikeind(peakind);

%now we'll step to either side of the peakspikeind and look for spikes that are
%separated by more than 1/threshold
if (opt.threshold > 0),
    dspike = 1/opt.threshold;
else
    %without a threshold, we'll just look for an empty bin and call that the end of the burst
    dspike = binsize;
end;
nspike = length(t);

nburst = length(peakspikeind);
on1 = zeros(1,nburst);
off1 = zeros(1,nburst);
for i = 1:nburst,
    a = peakspikeind(i);
    while ((a > 1) && (t(a)-t(a-1) <= dspike)),
        a = a-1;
    end;
    b = peakspikeind(i);
    while ((b < nspike) && (t(b+1)-t(b) <= dspike)),
        b = b+1;
    end;
    on1(i) = a;
    off1(i) = b;

    if (~isempty(opt.robustonoff) && (b-a > 2)),
        spikedt = median(diff(t(a:b)));
        if (t(a+1)-t(a) > opt.robustonoff*spikedt),
            on1(i) = a+1;
        end;
        if (t(b)-t(b-1) > opt.robustonoff*spikedt),
            off1(i) = b-1;
        end;
    end;
end;

%now check for the interburstdur
good = true(size(on1));
if (~isempty(opt.interburstdur) && (opt.interburstdur > 0)),
    tooclose = find(t(on1(2:end)) - t(off1(1:end-1)) <= opt.interburstdur);
    
    switch lower(opt.interburstaction)
      case 'merge',
        off1(tooclose) = off1(tooclose+1);
        good(tooclose+1) = false;
        
      case 'removeshortest',
        tooclosedur(1,:) = off1(tooclose) - on1(tooclose);
        tooclosedur(2,:) = off1(tooclose+1) - on1(tooclose+1);
        [q,sel] = min(tooclosedur);
        good(tooclose+sel-1) = false;
        
      case 'removeweakest',
        tooclosepk = [peakrate(tooclose); peakrate(tooclose+1)];
        [q,sel] = min(tooclosepk);
        good(tooclose+sel-1) = false;
    end;
end;
        
on1 = on1(good);
off1 = off1(good);
peakind = peakind(good);
peakrate = peakrate(good);

%merge bursts with overlapping onsets and offsets
if (opt.mergemultiples),
    overlap = off1(1:end-1) >= on1(2:end);
    good = [~overlap true];
    ismerged = [false overlap];
    
    ismerged = ismerged(good);
else
    %if we don't merge the multiples, then set the onsets and offsets so that they
    %don't overlap
    k = find(off1(1:end-1) >= on1(2:end));
    
    %for each overlap, look for the largest gap between spikes, and set the edges
    %around that gap
    for i = 1:length(k),
        b = peakspikeind(k(i)):peakspikeind(k(i)+1);
        ds = diff(t(b));
        
        [q,j] = max(ds);
        off1(k(i)) = b(j);
        on1(k(i)+1) = b(j)+1;
    end;
    good = true(size(on1));
    ismerged = false(size(on1));
end;

on = on1(good);
off = off1(good);
peakind = peakind(good);
peakrate = peakrate(good);

%figure out the burst centers - not necessarily the same as the index of the peak
burstctr = zeros(size(on));
for i = 1:length(on),
    k = on(i):off(i);
    burstctr(i) = mean(t(k));
end;

burstind = [on; off];




