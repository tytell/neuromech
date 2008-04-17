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

i = 1;
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
            
        otherwise,
            error('Unrecognized option %s', varargin{i});
    end;
end;

ispeak = false(size(y));
if (maxflat > 1),
    isflatpeak = false(size(y));
    flatpeakoffset = zeros(size(y));
end;
if (length(findpeaksign) > 1),
    peaksign = zeros(size(y));
end;

for i = 1:length(findpeaksign),
    ispeak1 = true(size(y));
    %can't detect peaks at the very beginning and end of the set
    ispeak1([1:numneighbors end-numneighbors+1:end]) = false;

    y1 = findpeaksign(i) * y;
    
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
        isflatpeak1 = true(size(y));
        isflatpeak1([1:k(1)-1 k(end)+1:end]) = false;
        
        isflatpeak1(k) = isflatpeak1(k) & (y1(k) == y1(k+flat));
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
peakind = find(ispeak);
if (maxflat > 1),
    isflatind = isflatpeak(ispeak);
    peakind(isflatind) = peakind(isflatind) + floor(flatpeakoffset(isflatpeak));
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
