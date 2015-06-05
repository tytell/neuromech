function burst = findbursts_smooth(t,amp, varargin)

opt.bindur = 0.02;

if (numel(varargin) == 1) && (isstruct(varargin{1}))
    opt = varargin{1};
else
    opt = parsevarargin(opt,varargin, 3);
    
    dt = t(2)-t(1);
    if (any(abs((diff(t)-dt)/dt) > 0.01)),
        opt.isSpikes = true;
    else
        opt.isSpikes = false;
    end;
end

if ~opt.isSpikes
    error('No supported yet');
end

tbin = floor(min(t)/opt.bindur)*opt.bindur : opt.bindur : (ceil(max(t)/opt.bindur)*opt.bindur);
[nspikes,binind] = histc(t,tbin);

spikerate = nspikes/opt.bindur;
spikeamp = accumarray(binind,abs(amp),size(spikerate),@nanmean,0);

[pkamp,pkind] = findpeaks2(amps);

for i = 1:length(pkind)
    a = pkind(i)-1;
    b = pkind(i)+1;
    while (a > 1) && (amp(a) > 0.5*pkamp(i))
        a = a-1;
    end
    while (b < length(amps)) && (amp(b) > 0.5*pkamp(i))
        b = b+1;
    end
    
end






