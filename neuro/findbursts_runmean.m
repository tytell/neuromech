function [burstctr burstdev burstskew spikectr tctr] = findbursts_runmean(t,spike,opt)

%if we're weighting by spike height, scale so that the weights are similar
%over the whole set
if (opt.isWeighted),
    spike = abs(spike);
    upper = 2.5*prctile(spike,75);
    spike(spike > upper) = upper;
end;

rundur2 = opt.rundur/2;

stepSize = opt.rundur*opt.stepFraction;        
tctr = (t(1)+rundur2):stepSize:(t(end)-rundur2);
if (~opt.isSpikes),
    dt = t(2)-t(1);
    tctrind = round((tctr-t(1))/dt) + 1;
    runlen2 = ceil(rundur2/dt);
else
    runind = [1 1];
end;

spikectr = repmat(NaN,size(tctr));
spikedev = repmat(NaN,size(tctr));
spikeskew = repmat(NaN,size(tctr));

for i = 1:length(tctr),
    if (opt.isSpikes),
        %find spikes within rundur2 of the current one
        while ((tctr(i) - t(runind(1))) > rundur2),
            runind(1) = runind(1)+1;
        end;
        while ((runind(2) <= length(t)) && (t(runind(2)) - tctr(i)) < rundur2),
            runind(2) = runind(2)+1;
        end;
        runind(2) = runind(2)-1;
    else
        runind = tctrind + [-runlen2 runlen2];
    end;
    
    %calculate the mean time, assuming there are any spikes in this
    %region
    %also a measure of uniformity and centrality - in a burst,
    %the spikes should be roughly centered in a single cluster
    %for the moment we'll use mean absolute deviation, weighted
    %by height
    %and skewness - log ratio of the mean spike deviation in the right
    %half of the region to the negative mean spike deviation in the
    %left half
    if (runind(1) <= runind(2)),
        k = runind(1):runind(2);

        dev = t(k) - tctr(i);
        right = dev > 0;
        if (opt.isWeighted),
            weight = sum(spike(k));
            spikectr(i) = sum(t(k).*spike(k)) ./ weight;
            spikedev(i) = sum(abs(dev).*spike(k)) ./ weight;
            
            if (sum(~right) <= 1)
                spikeskew(i) = Inf;
            elseif (sum(right) <= 1),
                spikeskew(i) = -Inf;
            else
                rightval = sum(dev(right).*spike(k(right))) ./ sum(spike(k(right)));
                leftval = -sum(dev(~right).*spike(k(~right))) ./ sum(spike(k(~right)));
                spikeskew(i) = log(rightval / leftval);
            end;
        else
            n = (runind(2)-runind(1)+1);
            spikectr(i) = sum(t(k)) ./ n;
            spikedev(i) = sum(abs(dev)) ./ n;
            if (all(right))
                spikeskew(i) = Inf;
            elseif (all(~right)),
                spikeskew(i) = -Inf;
            else
                rightval = sum(dev(right)) ./ sum(right);
                leftval = -sum(dev(~right)) ./ sum(~right);
                spikeskew(i) = log(rightval / leftval);
            end;
        end;
    end;
end;

%now look for places where the center of the spikes is relatively
%close to the center of the region over which the mean was taken.
%as the region passes over a burst, initially the burst center will
%be advanced relative to the region center, then it will
%progressively lag further and further, until it jumps to the next
%burst
dctr = spikectr - tctr;

%we want zero crossings that are decreasing.  Look at 4 points.
%Middle two should have opposite signs, and the ends should be
%further from zero than the middle
isBurst = false(size(dctr));
isBurst(2:end-2) = (sign(dctr(3:end-1)) <= 0) & (sign(dctr(2:end-2)) >= 0) & ...
    (dctr(1:end-3) > dctr(2:end-2)) & (dctr(3:end-1) > dctr(4:end));

%we want the point in which dctr is closest to zero (either
%positive or negative) as the best estimate of the burst center
burstind = find(isBurst);

[q,r] = min(abs([dctr(burstind); dctr(burstind+1)]));

%r is the row number of the minimum absolute value.  It'll be
%either 1 or 2.  By subtracting 1, we can use it as an offset to
%get the minimum dctr
burstctr = spikectr(burstind+r-1);
burstdev = spikedev(burstind+r-1);
burstskew = spikeskew(burstind+r-1);
