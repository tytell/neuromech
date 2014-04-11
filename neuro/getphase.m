function [phase,cycleind,cyclet] = getphase(t,sig, varargin)

opt.smoothparam = 1;    % lower is smoother
opt.threshold = [];     % peak threshold
opt = parsevarargin(opt,varargin,3);

if (~isempty(opt.threshold))
    thresh = {'threshold',opt.threshold};
else
    thresh = {};
end

%find the peaks, anticipating that some peaks may be flattened by clipping
nflat = 1 / (t(2) - t(1));      % 1 sec
[hival,hi] = findpeaks2(sig,'max', 'numneighbors',3, 'maxflat',nflat, thresh{:});
[loval,lo] = findpeaks2(sig,'min', 'numneighbors',3, 'maxflat',nflat, thresh{:});

npeak = min(length(hival),length(loval));
med = mean((hival(1:npeak) + loval(1:npeak))/2);

%also find the zeros
[zc,sgn] = findzeros(sig-med);

%get the ascending and descending zeros
asc = zc(sgn == 1);
desc = zc(sgn == -1);

%put them all together and sort them in order
pts = cat(1,asc(:), hi(:), desc(:), lo(:));
[pts,ord] = sort(pts(:));
%tp is the event type, 1=ascending zero, 2=maximum, 3=descending zero, 4=minimum
tp = catuneven(1,ones(length(asc),1), 2*ones(length(hi),1), ...
      3*ones(length(desc),1), 4*ones(length(lo),1));
tp = tp(ord);

%check to make sure we have at least two cycles to work with
if ((length(asc) < 2) || (length(desc) < 2) || (length(hi) < 2) || (length(lo) < 2))
    phase = NaN(size(t));
    cycleind = zeros(size(t));
    return;
end;
    
%iterate through the data set, trying to make it match the normal sequence:
%ascending zero, maximum, descending zero, minimum
good = true(size(pts));
ngood = sum(good);
nprevgood = 0;
iter = 1;

%first throw out repeated peaks or repeated zeros
while ((iter <= 10) && (ngood ~= nprevgood)),
    nprevgood = ngood;
    tp1 = tp(good);
    pts1 = pts(good);
    good1 = true(size(pts1));
    
    for i = 1:length(pts1)-1,
        if (~good1(i))
            continue;
        end;
        
        if (((tp1(i) == 2) || (tp1(i) == 4)) && ...
            ((tp1(i+1) == 2) || (tp1(i+1) == 4))),
            %keep the larger (absolute value) of two repeated peaks
            if (abs(sig(pts1(i))) > abs(sig(pts1(i+1))))
                good1(i+1) = false;
            else
                good1(i) = false;
            end;
        elseif (((tp1(i) == 1) || (tp1(i) == 3)) && ...
            ((tp1(i+1) == 1) || (tp1(i+1) == 3))),
            %keep the closer to zero of two repeated zeros
            if (abs(sig(pts1(i))) < abs(sig(pts1(i+1))))
                good1(i+1) = false;
            else
                good1(i) = false;
            end;
        end;
    end;
    good(good) = good1;
    ngood = sum(good);
    iter = iter + 1;
end;

pts = pts(good);
tp = tp(good);

%now run through the events and force them to match the correct sequence by
%adding in "missed" peaks or zeros as necessary
correctseq = tp(1);
k = 1;
pts2 = pts;
tp2 = tp;
for i = 1:length(tp),
    if (tp(i) ~= correctseq),
        %we're out of order, so insert the shortest sequence that will get us to the
        %right order
        if (tp(i) > tp(i-1))
            ins = (tp(i-1)+1:tp(i)-1)';
        else
            ins = [tp(i-1)+1:4 1:tp(i)-1]';
        end;
        
        tp2 = [tp2(1:k-1); ins; tp2(k:end)];
        pts2 = [pts2(1:k-1); zeros(size(ins)); pts2(k:end)];
        k = k+length(ins);
        correctseq = tp(i);
    end;
    
    %update the correct sequence value
    if (correctseq < 4)
        correctseq = correctseq+1;
    else
        correctseq = 1;
    end;
    k = k+1;
end;

tp = tp2;
pts = pts2;
good = pts ~= 0;

%estimate the number of cycles
dcycle = zeros(size(tp));
dcycle(tp == 1) = 1;
cyclenum = cumsum(dcycle);
if (cyclenum(1) == 0)
    cyclenum = cyclenum + 1;
end;
ncycle = cyclenum(end);

%get the time that corresponds to each event
cyclet = NaN(4,ncycle);
cycleind = zeros(4,ncycle);
ind = sub2ind(size(cyclet),tp(good),cyclenum(good));
cyclet(ind) = t(pts(good));
cycleind(ind) = pts(good);

ngoodcycle = sum(all(isfinite(cyclet)));
if (ngoodcycle > 1)
    %estimate the period, trying to correct for missed cycles as necessary
    dcycle = ones(1,ncycle);
    for iter = 1:2,
        cyclenum = cumsum(dcycle);
        
        per = NaN(size(cyclet));
        for i = 1:4,
            good = isfinite(cyclet(i,:));
            if (sum(good) > 1)
                per(i,good) = [diff(cyclet(i,good),[],2)./diff(cyclenum(good)) NaN];
            end;
        end;
        perall = nanmedian(per(:));
        
        %dcycle is at least one, but might be greater than one if we skip a cycle
        dcycle = round(nanmedian(per)./perall);
        dcycle(dcycle < 1) = 1;
    end;
    
    %estimate the period
    per = nanmedian(per,1);
    
    %estimate the phase for each event (ideally, it would be 0, 0.25, 0.5, 0.75, but the
    %waveform may not be so regular)
    ph1 = cyclet - cyclet(ones(4,1),:);
    ph1 = ph1 ./ per(ones(4,1),:);
    ph1 = nanmedian(ph1,2);
    ph1 = ph1(:,ones(1,ncycle));
    ph1 = ph1 + repmat(1:ncycle,[4 1]);
    
    good1 = isfinite(cyclet) & isfinite(ph1);
    good2 = (t >= min(cyclet(:,1))) & (t <= max(cyclet(:,end)));
    
    %interpolate the phase values
    phase = NaN(size(t));
    if (sum(good1(:)) > 3)
        h = nanmean(diff(cyclet(good1)));
        sp = csaps(cyclet(good1),ph1(good1), 1/(1+h^3/opt.smoothparam));
        phase(good2) = fnval(sp, t(good2));
    end
    
    %interpolate the zero crossings
    good = cycleind(1,:) ~= 0;
    toff = ((t(cycleind(1,good)+1) - t(cycleind(1,good))) ./ ...
        (sig(cycleind(1,good)+1) - sig(cycleind(1,good)))) .* ...
        (med - sig(cycleind(1,good)));
    cyclet(1,good) = cyclet(1,good) + toff;
    good = cycleind(3,:) ~= 0;
    toff = ((t(cycleind(3,good)+1) - t(cycleind(3,good))) ./ ...
        (sig(cycleind(3,good)+1) - sig(cycleind(3,good)))) .* ...
        (med - sig(cycleind(3,good)));
    cyclet(3,good) = cyclet(3,good) + toff;
else
    phase = NaN(size(t));
    cycleind = zeros(size(t));
end;


