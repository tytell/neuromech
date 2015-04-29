function burst = findbursts_surprise(t,spike, varargin)

opt.cycledur = 1;
opt.surprisethresh = 10;
opt.surprisedropthresh = 0.2;
opt.maxburstdur = 0.3;
opt.avgspikeratedur = 3;

if (numel(varargin) == 1) && (isstruct(varargin{1}))
    opt = varargin{1};
else
    opt = parsevarargin(opt,varargin, 3);
end

spikerate = zeros(size(t));
spikerate(2:end-1) = 1./((t(3:end) - t(1:end-2))/2);

avgspikerate0 = length(t) / (t(end) - t(1));

if ~isempty(opt.avgspikeratedur)
    avgspikeratedur2 = opt.avgspikeratedur/2;
end

%find the number of spikes in one second that would give us a surprise
%value above the threshold
Pthresh = exp(-opt.surprisethresh);
upperspikerate = poissinv(1-Pthresh,avgspikerate0);

%look for peaks in spike rate above that threshold
[~,pkind] = findpeaks2(spikerate, 'threshold',upperspikerate);

if isempty(pkind)
    warning('No bursts found');
end

k = 1;
for i = 1:length(pkind)
    a = pkind(i)-1;
    b = pkind(i)+1;
    
    if ~isempty(opt.avgspikeratedur)
        isavg = (t >= t(pkind(i))-avgspikeratedur2) & (t <= t(pkind(i))+avgspikeratedur2);
        tavg1 = t(find(isavg,1,'first'));
        tavg2 = t(find(isavg,1,'last'));
        avgspikerate1 = (sum(isavg)-1) / (tavg2 - tavg1);
    else
        avgspikerate1 = avgspikerate0;
    end
    
    s1 = surprisefcn(avgspikerate1,t(b)-t(a),2);
    sprev = s1;
    if (s1 >= opt.surprisethresh)
        good = true;
        while good && (t(b)-t(a) <= opt.maxburstdur)
            sleft = 0;
            sright = 0;
            sboth = 0;
            
            if (a > 1)
                sleft = surprisefcn(avgspikerate1,t(b)-t(a-1),b-a+1);
            end
            if (b < length(t))
                sright = surprisefcn(avgspikerate1,t(b+1)-t(a),b-a+1);
            end
            
            if (sleft >= opt.surprisethresh) && (sright >= opt.surprisethresh)
                sboth = surprisefcn(avgspikerate1,t(b+1)-t(a-1),b-a+2);
            end

            [snew,stepdir] = max([sboth sleft sright]);
            
            %if the new suprise is greater than the threshold and it
            %drops by less than the drop threshold, then we take it
            if (snew >= opt.surprisethresh) && ...
                    ((sprev-snew)/sprev < opt.surprisedropthresh)
                sprev = s1;
                s1 = snew;

                switch stepdir
                    case 1      % step both directions
                        a = a-1;
                        b = b+1;
                        
                    case 2      % step left
                        a = a-1;
                        
                    case 3      % step right
                        b = b+1;
                end
            else
                good = false;
            end
        end
        burst(k).surprise = s1;
        burst(k).avgspikerate = avgspikerate1;
        burst(k).spikerng = [a b];
        burst(k).trng = t([a b]);
        k = k+1;
    end
end

%merge overlapping bursts
burst0 = burst;

done = false;
i = 1;
while (i < 10) && ~done
    good = true(size(burst));
    for k = 1:length(burst)-1
        if good(k) && (burst(k).spikerng(2) >= burst(k+1).spikerng(1))
            if (all(burst(k).spikerng == burst(k+1).spikerng))
                %they're identical, so kill the second copy
                good(k+1) = false;
            elseif ((burst(k).spikerng(1) <= burst(k+1).spikerng(1)) && ...
                    (burst(k).spikerng(2) >= burst(k+1).spikerng(2)))
                %next burst is completely contained within the current one
                good(k+1) = false;
            elseif ((burst(k+1).spikerng(1) <= burst(k).spikerng(1)) && ...
                    (burst(k+1).spikerng(2) >= burst(k).spikerng(2)))
                %current burst is completely contained within the next one
                good(k) = false;
            else
                avgspikerate1 = 0.5*(burst(k).avgspikerate + burst(k+1).avgspikerate);
                a = burst(k).spikerng(1);
                b = burst(k+1).spikerng(2);
                
                if (t(b)-t(a) < opt.maxburstdur)
                    smerge = surprisefcn(avgspikerate1, t(b)-t(a), b-a);
                else
                    smerge = 0;
                end
                
                if (smerge >= opt.surprisethresh)
                    burst(k).surprise = smerge;
                    burst(k).avgspikerate = avgspikerate1;
                    burst(k).spikerng = [a b];
                    burst(k).trng = t([a b]);
                    
                    good(k+1) = false;
                else
                    dur1 = diff(burst(k).trng);
                    dur2 = diff(burst(k+1).trng);
                    s1 = burst(k).surprise;
                    s2 = burst(k).surprise;
                    
                    if (s1 >= s2) && (dur1 >= dur2)
                        good(k+1) = false;
                    elseif (s1 < s2) && (dur1 < dur2)
                        good(k) = false;
                    elseif (dur1 >= dur2) && ((s2-s1)/s2 < opt.surprisedropthresh)
                        good(k+1) = false;
                    elseif (dur2 > dur1) && ((s1-s2)/s1 < opt.surprisedropthresh)
                        good(k) = false;
                    elseif (s1 >= s2)
                        good(k+1) = false;
                    else
                        good(k) = false;
                    end
                end
            end
        end
    end
    done = all(good);
    burst = burst(good);
    i = i+1;
end

function S = surprisefcn(avgspikerate,T,nspikes)

nspikesexpected = avgspikerate*T;

P = poisscdf(nspikes,nspikesexpected,'upper');
S = -log(P);


        
