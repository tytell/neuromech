function [f,t] = filterspikes(b,a,spiket,sampfreq,method)
% function [f,t] = filterspikes(b,a,spiket,sampfreq,method)

opt.windowmin = 0.01;
opt.windowresamplefrac = 2;

%construct the time variable
t0 = min(spiket(:));
t1 = max(spiket(:));

dt = 1/sampfreq;
dthi = dt / opt.windowresamplefrac;

t = (0:dt:(t1-t0))';
thi = (0:dthi:(t1-t0))';

%check to see if we have channels in rows or columns
if (size(spiket,1) < size(spiket,2)),
    spiket = spiket';
    istranspose = true;
else
    istranspose = false;
end;

nchan = size(spiket,2);
nspike = size(spiket,1);

%turn b and a into cell arrays if they're not already
iscellout = iscell(b);
if (~iscell(b)),
    b = {b};
end;
if (~iscell(a)),
    a = {a};
    a = repmat(a,size(b));
end;
nwind = length(b);

spiket = spiket - t0;

f = cell(size(b));
switch lower(method),
    case 'bin',
        %bin the number of spikes to downsample, then run that data through
        %filtfilt
        binspike = histc(spiket,t);
    
        for i = 1:nwind,
            if (numel(b{i}) == 1),
                %build Gaussian windows
                windlen = b{i}(1) * sqrt(-log(opt.windowmin));
                
                twind = (-windlen:dt:windlen)';
                twind = twind - ((twind(1)+twind(end))/2);
                
                b1 = exp(-(twind/b{i}(1)).^2);
                b1 = b1 / sum(b1);
            end;

            f{i} = filtfilt(b1,a{i},binspike);
        end;
    
    case 'exact',
        if (numel(a{1}) ~= 1),
            error('Use of the a parameter is not implemented in exact mode');
        end;
        if (any(cellfun(@numel,b) ~= 1)),
            error('b parameter only specifies Gaussian window width.  Arbitrary filters not yet implemented.');
        end;
        
        for i = 1:nwind,
            windlen = b{i}(1) * sqrt(-log(opt.windowmin));
            twind = 0:dt:(windlen+dt);
            twind = [-twind(end:-1:2) twind];
            
            nwind = length(twind);
            iwind = 1:nwind;
            iwind = iwind - ((iwind(1)+iwind(end))/2);
            
            f1 = zeros(length(t),nchan);
            %get those spikes
            for ch = 1:nchan,                
                %exclude the spikes within half the window width at the
                %beginning and end of the time series
                good = isfinite(spiket(:,ch)) & ...
                    (spiket(:,ch) > -twind(1)+1/sampfreq) & ...
                    (spiket(:,ch) < t(end)-twind(end));
                
                spiket1 = spiket(good,ch);
                
                spikeind = floor(spiket1/dt) + 1;
                off = spiket1 - (spikeind-1)*dt;
                
                tt = repmat(twind,[length(spiket1) 1]) + repmat(off,[1 nwind]);
                ind = repmat(iwind,[length(spiket1) 1]) + repmat(spikeind,[1 nwind]);

                G = dt/(b{i}(1)*sqrt(pi))*exp(-(tt / b{i}(1)).^2);
                
                f1(:,ch) = accumarray(ind(:),G(:),[length(t) 1]);
            end;
            f{i} = f1;
        end;        
end;

t = t + t0;

if (istranspose),
    f = cellfun(@transpose,f,'UniformOutput',false);
    t = t';
end;
if (~iscellout),
    f = f{:};
end;





