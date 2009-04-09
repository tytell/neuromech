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


findpeaksign = 1;
numneighbors = 1;
maxflat = 1;
isstrict = false;
sortorder = 'none';
thresh = -Inf;
minpeakdistance = 0;
bias = 'morepeaks';

if (isnumeric(varargin{1}))
    dim = varargin{1};
    i = 2;
else
    dim = 1;
    i = 1;
end;

while (i <= length(varargin)),
    switch lower(varargin{i}),
      case 'max',
        findpeaksign = 1;
        i = i+1;
        
      case 'min',
        findpeaksign = -1;
        i = i+1;
        
      case 'minmax',
        findpeaksign = [1 -1];
        i = i+1;

      case 'thresh',
        thresh = varargin{i+1};
        i = i+2;
        
      case 'numneighbors',
        numneighbors = varargin{i+1};
        i = i+2;
        
      case 'maxflat',
        maxflat = varargin{i+1};
        i = i+2;
        
      case 'strict',
        isstrict = true;
        i = i+1;

      case 'sort',
        sortorder = varargin{i+1};
        i = i+2;

      case 'minpeakdistance',
        minpeakdistance = varargin{i+1};
        i = i+2;

      case {'morepeaks','biggerpeaks'},
        bias = lower(varargin{i});
        i = i+1;

      otherwise,
        error('Unrecognized option %s', varargin{i});
    end;
end;

ispeak = false(size(y));
isflatpeak = false(size(y));
flatpeakoffset = zeros(size(y));

if (length(findpeaksign) > 1),
    peaksign = zeros(size(y));
end;

for i = 1:length(findpeaksign),
    %adjust for minima or maxima
    thresh1 = findpeaksign(i) * thresh;
    y1 = findpeaksign(i) * y;

    %make sure we're over the threshold
    ispeak1 = y1 >= thresh1;
    %can't detect peaks at the very beginning and end of the set
    ispeak1([1:numneighbors end-numneighbors+1:end]) = false;

    %run through the neighbors
    k = numneighbors+1:length(y)-numneighbors;
    for off = 1:numneighbors,
        ispeak1(k) = ispeak1(k) & (y1(k) > y1(k+off)) & (y1(k) > y1(k-off));
        if (isstrict && (off > 1)),
            ispeak1(k) = ispeak1(k) & (y1(k+off-1) > y1(k+off)) & (y1(k-off+1) > y1(k-off));
        end;
    end;

    %run through the number of flat points
    for flat = 1:maxflat-1,
        k = numneighbors+1:length(y)-numneighbors-flat;
        isflatpeak1 = y1 >= thresh1;
        isflatpeak1([1:k(1)-1 k(end)+1:end]) = false;

        %to be a flat peak, all of the points at the top of the peak have to be equal
        for off = 1:flat,
            isflatpeak1(k) = isflatpeak1(k) & (y1(k) == y1(k+off));
        end;
        %and greater than the neighbors on either side
        for off = 1:numneighbors,
            isflatpeak1(k) = isflatpeak1(k) & (y1(k) > y1(k-off)) & (y1(k) > y1(k+flat+off));
            if (isstrict && (off > 1)),
                isflatpeak1(k) = isflatpeak1(k) & (y1(k-off+1) > y1(k-off)) & ...
                    (y1(k+flat+off-1) > y1(k+flat+off));
            end;
        end;
        
        isflatpeak = isflatpeak | isflatpeak1;
        flatpeakoffset(isflatpeak1) = flat/2;
        peaksign(isflatpeak1) = findpeaksign(i);
    end;

    if (length(findpeaksign) > 1),
        peaksign(ispeak1) = findpeaksign(i);
    end;
    ispeak = ispeak | ispeak1;
end;

if (maxflat > 1),
    ispeak = ispeak | isflatpeak;
end;
peaks = y(ispeak);
peakind = find(ispeak) + floor(flatpeakoffset(ispeak));

if (minpeakdistance > 1),
    %look for peaks that are separated by less than minpeakdistance, on either side
    good = true(size(peakind));
    good(1:end-1) = diff(peakind) >= minpeakdistance;

    %now find blocks of peaks, each of which is separated by less than minpeakdistance
    blockstart = makerow(find(good(1:end-1) & ~good(2:end))) + 1;
    blockend = makerow(find(~good(1:end-1) & good(2:end))) + 1;
    
    if (~isempty(blockstart) & ~isempty(blockend)),
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
        nmax = floor((peakind1(end)-peakind1(1))/minpeakdistance)+1;
        
        %sort the peaks from biggest to smallest
        p = peaks(block);
        [p,ord] = sort(p, 'descend');
        block = block(ord);
        peakind1 = peakind1(ord);
        
        good(block(1)) = true;        % always keep the biggest peak
        
        %step through and see if we can keep a few more peaks
        switch bias,
          case 'morepeaks',
            %bias towards more peaks means that we might throw out a larger peak for a
            %smaller one, if that means we can keep more
            keep = 1;
            b = 2;
            n = 1;
            while ((n < nmax) && (b <= length(block))),
                %if we find two that are separated by minpeakdistance, then keep them
                if (all(abs(peakind1(b) - peakind1(keep)) >= minpeakdistance)),
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
                   (peakind1(b) - peakind1(b-1) >= minpeakdistance)),
                good(block(b)) = true;
                n = n+1;
                b = b+1;
            end;
        end;
    end;
    
    peaks = peaks(good);
    peakind = peakind(good);
end;

switch lower(sortorder),
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
if ((length(findpeaksign) > 1) && (nargout == 3)),
    ps = peaksign(ispeak);
    if (~isempty(ord)),
        ps = ps(ord);
    end;
    varargout{3} = ps;
end;
