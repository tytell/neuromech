function varargout = analyzeKinematics(s,t,x,y, varargin)
% function [indpeak,confpeak,per,amp,midx,midy,exc,wavevel, ...
%			wavelen,waver,waven] = analyzeKinematics(s,t,x,y, good)
%	or		[per,amp,midx,midy,exc,wavevel,wavelen,waver,waven] = ...
%				analyzeKinematics(s,t,x,y, indpeak,[confpeak])
%	or		[indpeak,sgnpeak,per,amp,midx,midy,exc,wavevel, ...
%               wavelen,waver,waven] = ...
%				analyzeKinematics(s,t,x,y, indpeak,[confpeak], ...
%                                 'returnpeaksign')
%
% Pulls kinematic variables out of a midline.  Tracks body waves from
% tail to head and uses the positions to get the kinematic parameters.
% Fish must be swimming primarily along the x or y axes.
%
%	s - Distance along the midline.  Generally just linspace(0,fishlength,n),
%		where n is the number of points
%	t - Time, length = nfr
%	x,y - Midline position, size [n nfr].  (x(n,:),y(n,:)) is the tail
%		position
%	good - Points along the midline that you expect to get good body waves.
%		For n = 20, I usually use 10:20.
%
%	indpeak - Index of the time when the body wave peak is at a given point.
%		Has size [n npk], where npk is the number of peaks it managed to
%		track.  So indpeak(10,2) is the time index when peak 2 reached
%		point 10.  A good check after running the program is
%		plot(t,y(end,:), t(indpeak(end,:)),y(end,indpeak(end,:)),'o'),
%		which should put an o on the top of each peak in the tail beat
%	conf - No longer used, but still returned...
%	per, amp, wavelen - Period, amplitude, and wavelength, size [npt npk],
%		defined at each body point for each peak.
%	midx,midy - The "midline" of the midline.  Halfway between each side
%		to side peak.  Helps decide if the fish is swimming at an angle
%	exc - The excursion envelope.  Size [npt npk 2].  exc(:,:,1) is the
%		lateral excursion at each body point to one side and exc(:,:,2)
%		is the lateral excursion to the other side
%
% Note - you can edit indpeak, if it made errors, and rerun the program
% on it to calculate all the more useful parameters
%
% Mercurial revision hash: $Revision: 0570ba7d7b2b $ $Date: 2010/08/16 19:33:05 $
% Copyright (c) 2010, Eric Tytell

opt.allowreverse = true;
opt.nsmoothcurve = 5;
opt.dtsmoothcurve = [];
opt.dssmoothcurve = 0;
opt.curve = [];
opt.fwdwnd = 0.5;
opt.backwnd = 0.1;
opt.goodwavepts = [];
opt.peakmode = 'temporal';
opt.track = 'peaks';
opt.startpoint = 0.66;
opt.requiretail = false;
opt.returnpeaksign = false;

% Check parameters
if ((nargin >= 5) && (size(varargin{1},1) == length(s)) && (size(varargin{1},2) > 1)),
  isindpeak = true;
  indpeak = varargin{1};

  i = 2;
  good = true(size(indpeak,1),1);
elseif ((nargin >= 5) && ...
        (islogical(varargin{1}) || ...
         (isnumeric(varargin{1}) && ...
          all(size(varargin{1}) <= length(s)) && ...
          all((varargin{1} >= 1) & (varargin{1} <= length(s)))) ...
         )),
    good = varargin{1};
    i = 2;
    isindpeak = false;
else
    isindpeak = false;
    i = 1;
end;

if ((nargin == 4) || ~isnumeric(varargin{1}) ),
    good = true(size(s));
    good(1:4) = false;
    i = 1;
end;

opt = parsevarargin(opt,varargin,4+i, ...
                    'multival',{'track',{'peaks','zeros'}}, 'typecheck',false);

if (isnumeric(good)),
    k = good;
    good = false(size(s));
    good(k) = true;
end;
if (isempty(opt.goodwavepts)),
    opt.goodwavepts = true(size(s));
    opt.goodwavepts(3:end-2);
elseif (isnumeric(opt.goodwavepts)),
    g = false(size(s));
    g(opt.goodwavepts) = true;
    opt.goodwavepts = g;
end;

%estimate curvature
if (isempty(opt.curve)),
    curve = curvature(x,y,'spline','smooth',opt.dssmoothcurve);
else
    curve = opt.curve;
end;

%running mean
if (~isempty(opt.dtsmoothcurve)),
    opt.nsmoothcurve = floor(opt.dtsmoothcurve / (t(2)-t(1)));
end;
if (opt.nsmoothcurve > 0),
    curve = runavg(curve,opt.nsmoothcurve,2);
end;

nfr = size(x,2);
if (size(s,2) == nfr),
    error('s should be a vector');
end;
s = shiftdim(s);
npt = length(s);

if (~isindpeak),
    % Core of the procedure: Identify peaks/zeros in the curvature
    
    switch lower(opt.peakmode),
      case 'temporal',
        %look for peaks/zeros in time
        switch lower(opt.track),
          case 'peaks',
            ismax = false(size(curve));
            ismax(:,2:end-1) = (curve(:,2:end-1) >= curve(:,3:end)) & ...
                (curve(:,2:end-1) >= curve(:,1:end-2));
            ismin = false(size(curve));
            ismin(:,2:end-1) = (curve(:,2:end-1) <= curve(:,3:end)) & ...
                (curve(:,2:end-1) <= curve(:,1:end-2));

            ispeak = ismax | ismin;
            pksign = zeros(size(curve));
            pksign(ismax) = 1;
            pksign(ismin) = -1;
            
          case 'zeros',
            isasc = false(size(curve));
            % choose ascending points in which the second is closer to zero than the first
            isasc(:,2:end) = (curve(:,1:end-1) < 0) & (curve(:,2:end) > 0) & ...
                (curve(:,2:end) < -curve(:,1:end-1));
            % and ascending points in which the first is closer to zero than the second
            isasc(:,1:end-1) = isasc(:,1:end-1) | ...
                ((curve(:,1:end-1) < 0) & (curve(:,2:end) > 0) & ...
                (-curve(:,1:end-1) < curve(:,2:end)));

            isdesc = false(size(curve));
            isdesc(:,2:end) = (curve(:,1:end-1) > 0) & (curve(:,2:end) < 0) & ...
                (curve(:,2:end) < -curve(:,1:end-1));
            isdesc(:,1:end-1) = isdesc(:,1:end-1) | ...
                ((curve(:,1:end-1) > 0) & (curve(:,2:end) < 0) & ...
                (-curve(:,1:end-1) < curve(:,2:end)));

            ispeak = isasc | isdesc;
            pksign = zeros(size(curve));
            pksign(isasc) = 1;
            pksign(isdesc) = -1;
        end;
            
      case 'spatial',
        error('Spatial peak detection doesn''t work yet...');
        %look for peaks in space
        ismax = false(size(curve));
        ismax(2:end-1,:) = (curve(2:end-1,:) >= curve(3:end,:)) & ...
            (curve(2:end-1,:) >= curve(1:end-2,:));
        ismin = false(size(curve));
        ismin(2:end-1,:) = (curve(2:end-1,:) <= curve(3:end,:)) & ...
            (curve(2:end-1,:) <= curve(1:end-2,:));
        
      case 'both',
        error('Spatial/temporal peak detection doesn''t work yet...');
        ismax = false(size(curve));
        ismax(2:end-1,:) = (curve(2:end-1,:) >= curve(3:end,:)) & ...
            (curve(2:end-1,:) >= curve(1:end-2,:));
        ismax(:,2:end-1) = ismax(:,2:end-1) & (curve(:,2:end-1) >= curve(:,3:end)) & ...
            (curve(:,2:end-1) >= curve(:,1:end-2));
        ismin = false(size(curve));
        ismin(2:end-1,:) = (curve(2:end-1,:) <= curve(3:end,:)) & ...
            (curve(2:end-1,:) <= curve(1:end-2,:));
        ismin(:,2:end-1) = ismin(:,2:end-1) & (curve(:,2:end-1) <= curve(:,3:end)) & ...
            (curve(:,2:end-1) <= curve(:,1:end-2));
        
      otherwise,
        error('Unrecognized peakmode %s', opt.peakmode);
    end;
    
    %use the sign of the peak (ie, maximum vs minimum), not the sign of the curvature
    %itself.  That way we can track maxima, even if height of the maximum drops below
    %zero
    t0 = min(t) - 1;
    t = t - t0;		% make sure there's no t = 0

    if (strcmp(opt.track,'peaks') && ischar(opt.startpoint) && ...
        strcmp(opt.startpoint,'auto')),
        %find the point with the highest magnitude curvature, on average, that's also
        %symmetric left to right
        cpos = curve;
        cneg = curve;
        cpos(~ismax) = NaN;
        cneg(~ismin) = NaN;
        c = [nanmedian(cpos,2) -nanmedian(cneg,2)];
        
        rat = min(c,[],2) ./ max(c,[],2);
        choose = rat .* sum(c,2);
        
        [q,pt0] = max(choose);
    else
        pt0 = round(opt.startpoint * npt);
    end;
    
    %starting peaks are those we found for that point
    pk0 = find(ispeak(pt0,:));
    pksign0 = pksign(pt0,pk0);
    
    %search window
    wind = median(diff(pk0));
    wind = round(wind * [-opt.backwnd opt.fwdwnd]);
    wind = wind(1):wind(2);
    
    %check the starting peaks to make sure we don't have repeats
    for i = 1:length(pk0),
        d = pk0 - pk0(i);
        
        tooclose = (d >= wind(1)) & (d <= wind(end)) & ...
            (pksign0 == pksign0(i));
        if (sum(tooclose) > 1),
            %take the biggest one
            j = pk0(tooclose);
            [q,ind] = max(abs(curve(pt0,j)));
            ind = j(ind);
            
            pk0(tooclose) = NaN;
            pk0(i) = ind;
        end;
    end;

    pk0 = pk0(isfinite(pk0));
    
    npk = length(pk0);
    indpeak = zeros(npt,npk);
    indpeak(pt0,:) = pk0;
        
    %now step towards the tail, following each peak
    for i = 1:npk,
        curpk = pk0(i);
        curpkval = curve(pt0,pk0(i));
        curpksign = pksign(pt0,pk0(i));
        
        for pt = pt0+1:npt,
            wind1 = curpk + wind;
            wind1 = wind1((wind1 >= 1) & (wind1 <= nfr));
            
            %look for nearby peaks with the same sign curvature
            peaks = ispeak(pt,wind1) & (pksign(pt,wind1) == curpksign);
            peaks = wind1(peaks);
            
            if (~isempty(peaks)),
                %find the closest in time and height
                [q,ind] = min(abs(peaks - curpk) + abs(curve(pt,peaks) - curpkval));
                curpk = peaks(ind);
                curpkval = curve(pt,curpk);

                indpeak(pt,i) = curpk;
            else
                break;
            end;
        end;
    end;
    
    %now step towards the head
    for i = 1:npk,
        curpk = pk0(i);
        curpkval = curve(pt0,pk0(i));
        curpksign = pksign(pt0,pk0(i));
        
        for pt = pt0-1:-1:1,
            wind1 = curpk - wind;       %note the sign difference from above
            wind1 = wind1((wind1 >= 1) & (wind1 <= nfr));
            
            %look for nearby peaks with the same sign curvature
            peaks = find(ispeak(pt,wind1) & (pksign(pt,wind1) == curpksign));
            peaks = wind1(peaks);
            
            if (~isempty(peaks)),
                %find the closest in time and height
                [q,ind] = min(abs(peaks - curpk) + abs(curve(pt,peaks) - curpkval));
                curpk = peaks(ind);
                curpkval = curve(pt,curpk);

                indpeak(pt,i) = curpk;
            else
                break;
            end;
        end;
    end;

    indpeak(indpeak == 0) = NaN;

    % Make sure we tracked each peak for more than one time
    k = sum(isfinite(indpeak)) > 1;
    indpeak = indpeak(:,k);
end;

npk = size(indpeak,2);

goodpk = isfinite(indpeak);

% Get the actual times for each peak, rather than indices
tpeak = NaN(size(indpeak));
tpeak(goodpk) = t(indpeak(goodpk));

pt = repmat((1:npt)',[1 size(indpeak,2)]);
pkind = sub2ind(size(y),pt(goodpk),indpeak(goodpk));

% Get the positions for each peak
ypeak = NaN(size(indpeak));
ypeak(goodpk) = y(pkind);
xpeak = NaN(size(indpeak));
xpeak(goodpk) = x(pkind);
sgnpeak = zeros(size(indpeak));
sgnpeak(goodpk) = pksign(pkind);

% Sometimes we skip a peak.  If the sign of y at pt0 doesn't change
% from one peak to the next, then we skipped one.  Add columns of NaNs
% in between peaks where y has the same sign, to identify a missing peak
k = find((sgnpeak(pt0,1:end-1) ~= 0) & (sgnpeak(pt0,1:end-1) == sgnpeak(pt0,2:end)));
if (~isempty(k)),
    k = k + (1:length(k));
    npk2 = npk + length(k);
    k = setdiff(1:npk2,k);

    ip = NaN(npt,npk2);
    ip(:,k) = indpeak;
    indpeak = ip;

    yp = repmat(NaN,[npt npk2]);
    yp(:,k) = ypeak;
    ypeak = yp;
    xp = repmat(NaN,[npt npk2]);
    xp(:,k) = xpeak;
    xpeak = xp;
    tp = repmat(NaN,[npt npk2]);
    tp(:,k) = tpeak;
    tpeak = tp;
    
    npk = npk2;
end;

% Average number of frames it takes for a new peak to be formed at each body point
dp = nanmean(flatten(diff(indpeak,[],2)));

% Calculate wave length:
% Look for different peaks that are on the body at close to the same
% time
wavelen = NaN(size(indpeak));
for i = 1:size(indpeak,2)-1,
    % Time diference between this peak and the next one
    dist = repmat(indpeak(:,i)',[npt 1])-repmat(indpeak(:,i+1),[1 npt]);
    % Eliminate negative time differences, or ones that are greater
    % than dp
    dist((dist < -dp/2) | (dist > dp/2)) = NaN;

    % Get the minimum time difference.  These peaks are the ones that are
    % on the body at almost the same time
    [q,ind] = min(abs(dist));
    delta = dist(sub2ind(size(dist),ind,1:npt));
    
    % Linearly interpolate time to get the positions of the peaks at
    % the same time
    % but check to make sure the same peak wasn't found at the same frame time
    prevwavepos = s(ind);
    k = find(indpeak(ind+1,i+1) ~= indpeak(ind,i+1));
    %slope: change in s for a change in frame
    m = (s(ind(k)+1)-s(ind(k)))./(indpeak(ind(k)+1,i+1)-indpeak(ind(k),i+1));
    prevwavepos(k) = s(ind(k)) + delta(k)'.*m;

    % Wave length is twice the difference in peak postions at the same time
    wavelen(:,i) = 2*(s - prevwavepos);
end;

%check for empty wavelen values for the last tail beat
if (isnan(wavelen(end,end)) && isfinite(indpeak(end,end))),
    fr = indpeak(end,end);
    peaksign = sign(curve(end,fr));
    
    %guess the position of the previous peak
    guesspos = nanmedian(s(end) - (wavelen(:)/2)) ./ s(end) * length(s);
    
    cv1 = -peaksign * curve(:,fr);
    peakfr = find(cv1(2:end-1) > cv1(1:end-2) & ...
        cv1(2:end-1) > cv1(3:end)) + 1;
    
    [m,ind] = min(abs(peakfr - guesspos));
    ind = peakfr(ind);
    
    if (~isempty(ind)),
        wavelen(end,end) = 2*(s(end) - s(ind));
    else
        beep;
    end;
end;
   
% Calculate period: the time between every other peak
% Calculated in a complicated way because of skipped peaks
per(:,2:npk-1) = 2 * nanmean( cat(3, ...
    tpeak(:,3:end) - tpeak(:,2:end-1), ...
    tpeak(:,2:end-1) - tpeak(:,1:end-2)) ...
    ,3);
per(:,[1 npk]) = 2*(tpeak(:,[2 end]) - tpeak(:,[1 end-1]));

if (strcmp(opt.track,'peaks')),
    %Amplitude: half the mean distance between this peak and the previous one,
    %and between this peak and the next one
    amp(:,2:npk-1) = 0.5 * nanmean( cat(3, ...
                                        sqrt((ypeak(:,3:end)-ypeak(:,2:end-1)).^2 + (xpeak(:,3:end)-xpeak(:,2:end-1)).^2), ...
                                        sqrt((ypeak(:,2:end-1)-ypeak(:,1:end-2)).^2 + (xpeak(:,2:end-1)-xpeak(:,1:end-2)).^2)) ...
                                    ,3);
    amp(:,[1 npk]) = 0.5*sqrt((ypeak(:,[2 end])-ypeak(:,[1 end-1])).^2 + (xpeak(:,[2 end])-xpeak(:,[1 end-1])).^2);

    %Midway between each peak in x and y
    midy(:,2:npk) = 0.5*(ypeak(:,2:end) + ypeak(:,1:end-1));
    midy(:,1) = midy(:,2);
    midx(:,2:npk) = 0.5*(xpeak(:,2:end) + xpeak(:,1:end-1));
    midx(:,1) = midx(:,2);

    %Find missing columns and replace them with the next or the
    %previous value, depending on which is defined.
    [r,c] = find(isnan(midx));
    cplus = c+1;
    cminus = c-1;
    k = find((cminus >= 1) & (cplus <= npk));

    midx(sub2ind(size(midx),r(k),c(k))) = nanmean( cat(3, ...
                                                      midx(sub2ind(size(midx),r(k),cplus(k))), ...
                                                      midx(sub2ind(size(midx),r(k),cminus(k)))) ...
                                                   ,3);
    midy(sub2ind(size(midy),r(k),c(k))) = nanmean( cat(3, ...
                                                      midy(sub2ind(size(midy),r(k),cplus(k))), ...
                                                      midy(sub2ind(size(midy),r(k),cminus(k)))) ...
                                                   ,3);

    %Excursion: The distance between each peak position and the mid point
    exc = sqrt((ypeak - midy).^2 + (xpeak - midx).^2) .* sign(ypeak-midy);

    %Identify postive and negative peaks
    exchi = exc;
    exchi(exchi < 0) = NaN;
    exclo = exc;
    exclo(exclo > 0) = NaN;

    exc = cat(3,exclo,exchi);
else
    disp('amplitude measures don''t make much sense with zero tracking');
    amp = NaN(size(indpeak));
    midx = NaN(size(indpeak));
    midy = NaN(size(indpeak));
    exc = NaN([size(indpeak) 2]);
end;

% Calculate wave velocity: slope of the line through the peak position on
% the body plotted against the time it was at that position.  Only use
% body points identified as good
wavevel = NaN(1,npk);
waver = NaN(1,npk);
waven = NaN(1,npk);
for i = 1:size(tpeak,2),
    good1 = isfinite(tpeak(:,i)) & opt.goodwavepts;

    %If we have enough points, do the fit
    if (sum(good1) > 1),
        p = polyfit(tpeak(good1,i),s(good1),1);
        r = corrcoef(tpeak(good1,i),s(good1));

        wavevel(i) = p(1);
        waver(i) = r(1,2);
        waven(i) = sum(good1);
    else
        wavevel(i) = NaN;
        waver(i) = NaN;
        waven(i) = 0;
    end;
end;

% Remove columns with all NaNs or with NaNs for the tail
goodpk = any(isfinite(indpeak));
if (opt.requiretail)
    goodpk = goodpk & isfinite(indpeak(end,:));

    %except the last one
    l = last(any(isfinite(indpeak)) & ~isfinite(indpeak(end,:)));
    if (all(all(isnan(indpeak(:,l+1:end))))),
        goodpk(l) = true;
    end;
end;

indpeak = indpeak(:,goodpk);
conf = ones(size(indpeak));
per = per(:,goodpk);
amp = amp(:,goodpk);
midx = midx(:,goodpk);
midy = midy(:,goodpk);
exc = exc(:,goodpk,:);
wavelen = wavelen(:,goodpk);
wavevel = wavevel(goodpk);
waver = waver(goodpk);
waven = waven(goodpk);

if (isindpeak),
    varargout = {per,amp,midx,midy,exc,wavevel,wavelen,waver,waven};
else
    tpeak = tpeak + t0;
    if (opt.returnpeaksign),
        sgnpeak(sgnpeak == 0) = NaN;
        sgnpeak = mode(sgnpeak);
        varargout = {indpeak,sgnpeak, per,amp,midx,midy,exc,wavevel,wavelen,waver,waven};
    else
        varargout = {indpeak,conf, per,amp,midx,midy,exc,wavevel,wavelen,waver,waven};
    end;
end;

