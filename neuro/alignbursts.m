function [burstt,ord] = alignbursts(burstt0, varargin)
% function [burstt,ord] = alignBursts(burstt0, varargin)
% Align bursts from several channels, assuming that the bursts are coming
% at approximately the same frequency
% ord is the order in the new variable.  So if you have another variable
% like burstt0 (say burstper0), you can reorder it as follows:
%   burstper = nans(size(burstt));
%   burstper(ord) = cat(2,burstper0{:});

maxdiff = Inf;
channelorder = [];
nperavg = 4;
aligncycles = false;

i = 1;
while (i <= length(varargin)),
    switch lower(varargin{i}),
        case 'maxdiff',
            maxdiff = varargin{i+1};
            i = i+2;
        case 'channelorder',
            channelorder = varargin{i+1};
            i = i+2;
        case 'aligncycles',
            aligncycles = true;
            i = i+1;
            
        otherwise,
            error('Unrecognized option %s\n', varargin{i});
    end;
end;

nchan = length(burstt0);

chan = cell(size(burstt0));
for i = 1:nchan,
    chan{i} = i*ones(size(burstt0{i}));
end;

burstt1 = cat(2,burstt0{:});
chan = cat(2,chan{:});

%sort the burst times to match up corresponding bursts, and leave NaNs
%for skipped bursts
[burstt1,ord] = sort(burstt1);
chan = chan(ord);

if (isempty(channelorder)),
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
end;

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
burstt = nans(nchan,nbursts);
ord2 = sub2ind(size(burstt),chan,cols);
burstt(ord2) = burstt1;

%this allows us to do the sorting in time order (which gives us the ord
%matrix) and the sorting by channels (which gives us the ord2 matrix) all
%in one operation
[q,revord] = sort(ord);
ord = ord2(revord);

if (aligncycles),
    good = false(size(burstt));
    good(:,1:nperavg) = true;

    i = nperavg+1;
    while (i <= nbursts),
        meanburst = nanmean(burstt(:,i-nperavg:i));
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

        
    
