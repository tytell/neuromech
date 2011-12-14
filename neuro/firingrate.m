function varargout = firingrate(spiket, binsize, smooth, varargin)
% function rate = firingrate(spiket, binsize, smooth, method)
% smooth time interval over which to perform a moving average
% if binsize is zero and smooth is zero, returns the spike-wise instantaneous firing rate.
% if binsize is zero, but smooth is non-zero, returns the spike-wise instantaneous firing
% rate, smoothed using csaps and interpolated on to evenly spaced samples with spacing
% specified by resamp

opt.truncategaussian = 0.005;
opt.method = 'hist';

if (nargin > 3)
    opt = parsevarargin(opt,varargin,4, ...
        'multival',{'method',{'hist','spline','gaussian'}});
elseif (nargin == 2)
    smooth = 0;
end;

if (iscell(spiket)),
    for ch = 1:length(spiket),
        spiket{ch} = makecol(spiket{ch});
    end;
    
    spiket = catuneven(2,spiket{:});
end;

if (size(spiket,1) < size(spiket,2)),
    spiket = spiket';
    istransposed = true;
else
    istransposed = false;
end;

nchan = size(spiket,2);
iseven = true;
edges = [];

if (numel(smooth) == 1),
    smooth = smooth(ones(1,nchan),1);
end;

if (numel(binsize) > 1),
    edges = binsize(:);

    binsize = edges(2)-edges(1);
    if (all(abs((diff(edges)-binsize)/binsize) < 0.001)),
        iseven = true;
        binsize = edges(2)-edges(1);
    else
        iseven = false;
        binsize = NaN;
    end;
end;

switch lower(opt.method),
    case 'gaussian',
        rate0 = zeros(size(spiket));
        spikedt = diff(spiket);
        
        rate0(2:end-1,:) = 1 ./ ((spikedt(1:end-1) + spikedt(2:end))/2);
        t0 = spiket;
        
        if (isempty(edges))
            tstart = floor(min(spiket(:))/binsize)*binsize - binsize/2;
            tend = ceil(max(spiket(:))/binsize)*binsize + binsize/2;
            
            t = (tstart+binsize/2:binsize:tend)';
        else
            t = edges;
            tstart = edges(1);
        end;
        nbin = length(t);
        
        bin = floor((t0 - tstart)/binsize)+1;
        maxbinwidth = ceil(smooth * sqrt(-2*log(opt.truncategaussian)));
        
        rate = NaN(nbin,nchan);
        spikeindbybin = zeros(nbin,nchan);
        for c = 1:nchan,
            spikeindbybin1 = accumarray(bin(:,c),(1:size(spiket,1))',[nbin 1],@(x) {x(:)});
            
            for i = 1:nbin,
                k = (i-maxbinwidth):(i+maxbinwidth);
                k = k((k >= 1) & (k <= nbin));
                
                spikeind1 = cat(1,spikeindbybin1{k});
                
                if (~isempty(spikeind1))
                    spiket1 = spiket(spikeind1);
                    rate1 = rate0(spikeind1);
                    
                    weight = exp(-0.5*((spiket1 - t(i))/smooth).^2);
                    [~,j] = max(weight);
                    
                    rate(i,c) = sum(weight.*rate1)/sum(weight);
                    spikeindbybin(i,c) = spikeind1(j);
                else
                    rate(i,c) = 0;
                end;
            end;
            
        end;
        
        if (istransposed)
            spikeindbybin = spikeindbybin';
        end;
    
    case 'spline',
        rate0 = NaN(size(spiket));
        rate0(2:end-1,:) = 1 ./ ((spiket(3:end,:) - spiket(1:end-2,:))/2);
        t0 = spiket;
        
        if (nargin >= 3),
            t = t0;
            if (nargin == 4),
                tstart = floor(min(spiket(:))/binsize)*binsize - binsize/2;
                tend = ceil(max(spiket(:))/binsize)*binsize + binsize/2;
                
                t = (tstart:binsize:tend)';
            end;
            
            rate = NaN(length(t),size(rate0,2));
            for i = 1:size(rate,2),
                good = isfinite(rate0(:,i));
                a = ceil((first(t0(:,i),good) - tstart) / binsize) + 1;
                b = floor((last(t0(:,i),good) - tstart) / binsize);
                
                rate1 = windavg(t0(good,i),rate0(good,i), smooth(i)*binsize,t(a:b));
                rate1(isnan(rate1)) = 0;
                rate(a:b,i) = rate1;
            end;
        else
            t = t0;
            rate = rate0;
        end;
        
        bin = [];
        
    case 'hist',
        if (isempty(edges)),
            t0 = floor(min(spiket(:))/binsize)*binsize - binsize/2;
            t1 = ceil(max(spiket(:))/binsize)*binsize + binsize/2;
            
            edges = (t0:binsize:t1)';
            edges(end) = Inf;
        end;
        
        n = zeros(length(edges),nchan);
        
        bin = NaN(size(spiket));
        for ch = 1:nchan,
            good = isfinite(spiket(:,ch));
            [n(:,ch),bin1] = histc(spiket(good,ch),edges);
            bin(good,ch) = bin1;
        end;
        n = n(1:end-1,:);
        
        if (iseven),
            rate = n / binsize;
            t = edges(1:end-1) + binsize/2;
        else
            rate = n ./ diff(edges);
            t = (edges(1:end-1) + edges(2:end))/2;
        end;
        
        if (smooth ~= 0),
            for i = 1:nchan,
                rate(:,i) = runavg(rate(:,i),smooth(i));
            end;
        end;
end;

if (istransposed),
    t = t';
    rate = rate';
    bin = bin';
end;

switch nargout,
    case 1,
        varargout = {rate};
    case 2,
        varargout = {t,rate};
    case 3,
        varargout = {t,rate,bin};
    case 4,
        varargout = {t,rate,bin,spikeindbybin};
end;

    