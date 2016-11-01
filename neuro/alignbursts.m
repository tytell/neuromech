function [burstt,ord] = alignbursts(burstt0, varargin)
% function [burstt,ord] = alignBursts(burstt0, varargin)
% Align bursts from several channels, assuming that the bursts are coming
% at approximately the same frequency
% ord is the order in the new variable.  So if you have another variable
% like burstt0 (say burstper0), you can reorder it as follows:
%   burstper = nans(size(burstt));
%   burstper(ord) = cat(2,burstper0{:});

opt.maxdiff = Inf;
opt.channelorder = [];
opt.nperavg = 4;
opt.aligncycles = false;

opt = parsevarargin(opt,varargin, 2);


if iscell(burstt0)
    nchan = length(burstt0);
    
    chan = cell(size(burstt0));
    for i = 1:nchan,
        chan{i} = i*ones(size(burstt0{i}));
    end;

    burstt1 = cat(2,burstt0{:});
    chan = cat(2,chan{:});
else
    nchan = size(burstt0,2);
    
    chan = repmat(1:nchan,[size(burstt0,1) 1]);
    
    burstt1 = burstt0(:)';
    chan = chan(:)';
end

%sort the burst times to match up corresponding bursts, and leave NaNs
%for skipped bursts
[burstt1,ord] = sort(burstt1);
chan = chan(ord);

if (isempty(opt.channelorder)),
    ischan = false(nchan,1);
    channelorder = zeros(nchan,1);
    j = 1;
    i = 1;
    while (i <= length(chan)) && (j < nchan),
        if (~ischan(chan(i))),
            channelorder(j) = chan(i);
            j = j+1;
            ischan(chan(i)) = true;
        end;
        i = i+1;
    end;
    lastchan = find(~ischan);
    channelorder(j:end) = lastchan;
else
    channelorder = opt.channelorder;
end

%now go through the sorted channel orders and look for any instances in
%which the channels repeat in the wrong order - that means we got a double 
%burst on that channel or a skipped burst on another
cols = zeros(size(chan));       % column numbers
chanstep = 0;
col = 1;
for i = 1:length(chan),
    nextstep = find(channelorder == chan(i));
    if (nextstep <= chanstep),
        col = col+1;
    end;

    cols(i) = col;
            
    chanstep = nextstep;
end;

%sort the channels so that each column contains bursts that match up
%according to the procedure above
nbursts = max(cols);
burstt = NaN(nchan,nbursts);
ord2 = sub2ind(size(burstt),chan,cols);
burstt(ord2) = burstt1;

%this allows us to do the sorting in time order (which gives us the ord
%matrix) and the sorting by channels (which gives us the ord2 matrix) all
%in one operation
[q,revord] = sort(ord);
ord = ord2(revord);

if (opt.aligncycles),
    good = false(size(burstt));
    good(:,1:opt.nperavg) = true;

    i = opt.nperavg+1;
    while (i <= nbursts),
        meanburst = nanmean(burstt(:,i-opt.nperavg:i));
        per = diff(meanburst);
        meanper = nanmean(per);

        for j = 1:nchan,
            off = 1;
            err = Inf;
            done = false;
            while ((i+off <= nbursts) && ~done),
                curper = burstt(j,i+off) - burstt(j,i);
                if (isfinite(curper)),
                    if (abs(curper - meanper) > err),
                        done = true;
                    else
                        err = abs(curper - meanper);
                        off = off+1;
                    end;
                else
                    off = off+1;
                end;
            end;

            if (done),
                good(j,i+off) = true;
            end;
        end;
    end;
end;

        
    
