function varargout = findpeaks2(y, varargin)
% function [peak,peakind,...] = findpeaks2(y, ...)
% Finds peaks in y.
% Options:
%     'max','min',or 'minmax' - Finds maxima, minima or both (default = max)
%     'numneighbors' - Makes sure each peak is greater (less) than
%         numneighbors on both sides (default = 1)
%     'maxflat' - Maximum number of points at the top of a peak that can be
%         exactly equal to each other (default = 1)
%     'strict' - Strict peak finding, meaning that that the top of the peak
%         must be greater (less) than each neighbor, but also that each
%         neighbor must be greater (less) than the successive neighbor on
%         the peak's slope
%     'sort' - Sorts the peaks.  Options are 'up','down','absup','absdown'
%         or 'none'
%     'thresh' - Threshold for the peaks.  NB: for minmax, this is an absolute
%         threshold, but for the others, its sign is not changed
%     'minpeakdistance' - Minimum distance between peaks
%     'morepeaks','biggerpeaks' - Defines the bias for how to sort peaks that are too
%        close.  'morepeaks' will potentially throw out a large peak, if it allows the
%        minpeakdistance criterion to hold, while 'biggerpeaks' always takes the largest
%        peaks.  Default is 'morepeaks'
%
% Returns:
%     peak - The values of the peaks, in the same order as y.
%     peakind - The indices of the peaks.  For flat peaks, returns the
%        index of the central point (rounded down)
%     peaksign - Only returned for 'minmax' finding.  Returns 1 for maxima
%        and -1 for minima.
%
% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2010, Eric Tytell <tytell at jhu dot edu>


opt.findpeaks = 'max';
opt.findpeaksign = 1;
opt.numneighbors = 1;
opt.maxflat = 1;
opt.strict = false;
opt.sortorder = 'none';
opt.threshold = -Inf;
opt.minpeakdistance = 0;
opt.bias = 'morepeaks';

opt = parsevarargin(opt,varargin, 'multival', ...
    {'findpeaks',{'min','max','minmax'}, 'bias',{'morepeaks','biggerpeaks'}}, ...
    'synonyms',{{'sortorder',{'sort'}}});

switch opt.findpeaks,
    case 'max',
        opt.findpeaksign = 1;
    case 'min',
        opt.findpeaksign = -1;
    case 'minmax',
        opt.findpeaksign = [-1 1];
end;

ispeak = false(size(y));
isflatpeak = false(size(y));
flatpeakoffset = zeros(size(y));

if (length(opt.findpeaksign) > 1),
    peaksign = zeros(size(y));
end;

for i = 1:length(opt.findpeaksign),
    %adjust for minima or maxima
    thresh1 = opt.threshold;
    y1 = opt.findpeaksign(i) * y;

    %make sure we're over the threshold
    ispeak1 = y1 >= thresh1;
    %can't detect peaks at the very beginning and end of the set
    ispeak1([1:opt.numneighbors end-opt.numneighbors+1:end]) = false;

    %run through the neighbors
    k = opt.numneighbors+1:length(y)-opt.numneighbors;
    for off = 1:opt.numneighbors,
        ispeak1(k) = ispeak1(k) & (y1(k) > y1(k+off)) & (y1(k) > y1(k-off));
        if (opt.strict && (off > 1)),
            ispeak1(k) = ispeak1(k) & (y1(k+off-1) > y1(k+off)) & (y1(k-off+1) > y1(k-off));
        end;
    end;

    %run through the number of flat points
    flatdone = false;
    flat = 1;
    while (~flatdone && (flat <= opt.maxflat-1)),
        k = opt.numneighbors+1:length(y)-opt.numneighbors-flat;
        isflatpeak1 = y1 >= thresh1;
        isflatpeak1([1:k(1)-1 k(end)+1:end]) = false;

        %to be a flat peak, all of the points at the top of the peak have to be equal
        for off = 1:flat,
            isflatpeak1(k) = isflatpeak1(k) & (y1(k) == y1(k+off));
        end;
        
        %we can declare ourselves done if we don't have any runs of length flat
        flatdone = all(~isflatpeak1);
        
        %and greater than the neighbors on either side
        for off = 1:opt.numneighbors,
            isflatpeak1(k) = isflatpeak1(k) & (y1(k) > y1(k-off)) & (y1(k) > y1(k+flat+off));
            if (opt.strict && (off > 1)),
                isflatpeak1(k) = isflatpeak1(k) & (y1(k-off+1) > y1(k-off)) & ...
                    (y1(k+flat+off-1) > y1(k+flat+off));
            end;
        end;
        
        isflatpeak = isflatpeak | isflatpeak1;
        flatpeakoffset(isflatpeak1) = flat/2;
        peaksign(isflatpeak1) = opt.findpeaksign(i);
        
        flat = flat+1;
    end;

    if (length(opt.findpeaksign) > 1),
        peaksign(ispeak1) = opt.findpeaksign(i);
    end;
    ispeak = ispeak | ispeak1;
end;

if (opt.maxflat > 1),
    ispeak = ispeak | isflatpeak;
end;
peaks = y(ispeak);
peakind = find(ispeak) + floor(flatpeakoffset(ispeak));

if (opt.minpeakdistance > 1),
    %look for peaks that are separated by less than minpeakdistance, on either side
    good = true(size(peakind));
    good(1:end-1) = diff(peakind) >= opt.minpeakdistance;

    %now find blocks of peaks, each of which is separated by less than minpeakdistance
    blockstart = makerow(find(good(1:end-1) & ~good(2:end))) + 1;
    blockend = makerow(find(~good(1:end-1) & good(2:end))) + 1;
    
    if (~isempty(blockstart) && ~isempty(blockend)),
        if (blockend(1) < blockstart(1)),
            blockstart = [1 blockstart];
        end;
        if (blockstart(end) > blockend(end)),
            blockend = [blockend length(peakind)];
        end;
    end;
    
    %run through the blocks
    for i = 1:length(blockstart),
        block = blockstart(i):blockend(i);
        peakind1 = peakind(block);
        
        %initially, assume we're going to get rid of all of the peaks
        good(block) = false;

        %this is the largest number of peaks that can (potentially) fit in the space
        nmax = floor((peakind1(end)-peakind1(1))/opt.minpeakdistance)+1;
        
        %sort the peaks from biggest to smallest
        p = peaks(block);
        [p,ord] = sort(p, 'descend');
        block = block(ord);
        peakind1 = peakind1(ord);
        
        good(block(1)) = true;        % always keep the biggest peak
        
        %step through and see if we can keep a few more peaks
        switch opt.bias,
          case 'morepeaks',
            %bias towards more peaks means that we might throw out a larger peak for a
            %smaller one, if that means we can keep more
            keep = 1;
            b = 2;
            n = 1;
            while ((n < nmax) && (b <= length(block))),
                %if we find two that are separated by minpeakdistance, then keep them
                if (all(abs(peakind1(b) - peakind1(keep)) >= opt.minpeakdistance)),
                    good(block(b)) = true;
                    keep = [keep b];
                    b = b+1;
                    n = n+1;
                else
                    %otherwise, look at the next smaller peak
                    b = b+1;
                end;
            end;

          case 'biggerpeaks',
            %bias towards bigger peaks: step down the list of peaks in order of
            %height.  If they're separated by >= minpeakdistance, then keep them.  Stop
            %at the first pair separated by < minpeakdistance
            b = 2;
            n = 1;
            while ((n < nmax) && (b <= length(block)) && ...
                   (peakind1(b) - peakind1(b-1) >= opt.minpeakdistance)),
                good(block(b)) = true;
                n = n+1;
                b = b+1;
            end;
        end;
    end;
    
    peaks = peaks(good);
    peakind = peakind(good);
end;

switch lower(opt.sortorder),
    case {'none',''},
        ord = [];  
        
    case {'up','ascending'},
        [peaks,ord] = sort(peaks);
        peakind = peakind(ord);
        
    case {'down','descending'},
        [peaks,ord] = sort(-peaks);
        peaks = -peaks;
        peakind = peakind(ord);
        
    case {'absdown','absdescending'},
        [q,ord] = sort(-abs(peaks));
        peaks = peaks(ord);
        peakind = peakind(ord);
        
    case {'absup','absascending'},
        [q,ord] = sort(abs(peaks));
        peaks = peaks(ord);
        peakind = peakind(ord);
        
    otherwise,
        error('Unrecognized sort order');
end;

varargout = {peaks, peakind};
if ((length(opt.findpeaksign) > 1) && (nargout == 3)),
    ps = peaksign(ispeak);
    if (~isempty(ord)),
        ps = ps(ord);
    end;
    varargout{3} = ps;
end;
