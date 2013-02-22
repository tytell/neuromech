function varargout = eststimfreq(stimt,stim,protocol,varargin)

opt.method = 'peaks&zeros';
opt.smoothphasedur = 0.1;
opt.quiet = true;
opt.smoothparam = 1;            % lower is smoother

opt = parsevarargin(opt,varargin,2);

n = round(opt.smoothphasedur / (stimt(2) - stimt(1)));

if ((nargin == 2) || isempty(protocol))
    protocol.t = 0;
    protocol.duration = Inf;
elseif (isfield(protocol,'StartSamps'))
    protocol.t = double(protocol.StartSamps)/protocol.SampleRate;
    tend = double(protocol.EndSamps)/protocol.SampleRate;
    protocol.duration = tend - protocol.t;
end;

istrans = false;
if ((size(stimt,2) == 1) && (size(stim,2) == 1))
    istrans = true;
    stimt = stimt';
    stim = stim';
end

done = false;
if (isfield(protocol,'OutputType'))
    phase = NaN(size(stimt));
    freq = NaN(size(stimt));
    amp = NaN(size(stimt));

    if (iscellstr(protocol.OutputType) && (numel(protocol.OutputType) == 1))
        protocol.OutputType = protocol.OutputType{1};
    end;
    switch protocol.OutputType,
        case 'Sine',
            for i = 1:length(protocol.t),
                istreat = (stimt >= protocol.t(i)) & (stimt < protocol.t(i)+protocol.duration(i));
                
                freq(istreat) = protocol.Frequency(i)*protocol.SampleRate / (2*pi);
                phase(istreat) = (stimt(istreat)-protocol.t(i)).*freq(istreat);
                
%                 k = first(istreat)+(0:sum(istreat)-1);
%                 amp1 = interp1(phase(istreat), stim(istreat), 0:0.25:max(phase(istreat)));
%                 ampind = interp1(phase(istreat), k, 0:0.25:max(phase(istreat)));
%                 len = floor((length(amp1)-1)/4)*4 + 1;
%                 amp1 = amp1(1:len);
%                 ampind = ampind(1:len);
%                 
%                 ampval(3:4:end) = abs(amp1(2:4:end) - amp1(4:4:end));
%                 ampval(5:4:end) = abs(amp1(6:4:end) - amp1(4:4:end));
%                 
%                 amp(istreat) = interp1(ampind(3:2:end),ampval(3:2:end), k);
            end;
            done = true;
            
        case 'Frequency sweep',
            istreat = (stimt >= protocol.t) & (stimt < protocol.t+protocol.duration);
            
            k = protocol.SweepRate*protocol.SampleRate^2;
            f0 = protocol.Frequency*protocol.SampleRate;
            tt = stimt(istreat) - protocol.t;
            
            freq(istreat) = (tt*k + f0) / (2*pi);
            phase(istreat) = (0.5*k*tt.^2 + f0 * tt) / (2*pi);
            
            done = true;
    end;
elseif (isfield(protocol,'hzpersec'))
    phase = NaN(size(stimt));
    freq = NaN(size(stimt));
    amp = NaN(size(stimt));

    istreat = (stimt >= protocol.t(1)) & (stimt < protocol.t(1)+protocol.duration(1));
    
    k = protocol.hzpersec;
    f0 = protocol.frequency;
    tt = stimt(istreat) - protocol.t(1);
    
    freq(istreat) = (tt*k + f0);
    phase(istreat) = (0.5*k*tt.^2 + f0 * tt);
    
    done = true;
end;    

if (done)
    if (istrans)
        phase = phase';
        freq = freq';
        amp = amp';
    end
    
    varargout = {phase,freq,amp};
    return;
end;
            
switch opt.method,
  case 'phaser',
    phase = NaN(size(stimt));
    phase0 = NaN(size(stimt));

    for i = 1:length(protocol.t),
        if (~opt.quiet)
            fprintf('Treatment %d...\n',i);
        end;
        istreat = (stimt >= protocol.t(i)) & (stimt < protocol.t(i)+protocol.duration(i));
        
        stim1 = stim(istreat);
        if (~isempty(stim1)),
            [phr,phase1] = newPhaser(stim1);
            phase1 = phase1(:)' / (2*pi);

            if (opt.smoothphasedur > 0),
                phase0(istreat) = phase1;
                
                phase1 = runavg(phase1,n,2);
                phase(istreat) = phase1;
            else
                phase(istreat) = phase1;
            end;
        end;
    end;
    
    freq = deriv(stimt,phase);

    if (istrans)
        phase = phase';
        freq = freq';
        phase0 = phase0';
    end
    
    varargout = {phase,freq, phase0};
    
  case 'peaks&zeros',
    phase = NaN(1,length(stimt));
    cycleevent = zeros(1,length(stimt));
    amplitude = NaN(1,length(stimt));
    
    for i = 1:length(protocol.t),
        if (~opt.quiet)
            fprintf('Treatment %d...\n',i);
        end;
        istreat = (stimt >= protocol.t(i)) & (stimt < protocol.t(i)+protocol.duration(i));
        
        med = nanmedian(stim(istreat));
        stim1 = stim(istreat);
        if (~isempty(stim1)),
            [phase1,cycleind1] = peaksandzeros(stimt(istreat),stim1-med, opt.smoothparam);

            if (size(cycleind1,1) ~= 4)
                warning('Non-sinusoidal signal at treatment %d\n', i);
                continue;
            end;

            phase(istreat) = phase1;

            %estimate the stimulus amplitude
            A1 = NaN(size(cycleind1));
            good = cycleind1 ~= 0;
            A1(good) = stim1(cycleind1(good));
            
            amp1 = NaN(2,size(cycleind1,2));
            
            %amplitude is half the distance from a max to a min (or vice versa)
            amp1(1,:) = (A1(2,:) - A1(4,:)) / 2;
            amp1(2,1:end-1) = (A1(2,2:end) - A1(4,1:end-1))/2;

            %and it's defined at the point halfway in between the max and min
            ampind = zeros(size(amp1));
            ampind(1,:) = round((cycleind1(2,:) + cycleind1(4,:))/2);
            ampind(1,(cycleind1(2,:) == 0) | (cycleind1(4,:) == 0)) = 0;
            
            ampind(2,1:end-1) = round((cycleind1(4,1:end-1) + cycleind1(2,2:end))/2);
            ampind(2,[cycleind1(4,1:end-1)==0 false]) = 0;
            ampind(2,[cycleind1(2,2:end)==0 false]) = 0;
            
            %good points
            good1 = (ampind ~= 0) & ~isnan(amp1);
            a = min(ampind(ampind ~= 0));
            b = max(ampind(ampind ~= 0));

            %spline the amplitude
            amp2 = NaN(size(stim1));

            %but use a smoothing spline so we don't end up with weird overshoots in the
            %gaps
            h = nanmean(diff(ampind(good1)));
            sp = csaps(ampind(good1),amp1(good1), 1/(1+h^3/opt.smoothparam));
            
            amp2(a:b) = fnval(sp, a:b);
            amplitude(istreat) = amp2;
            
            a = first(istreat);
            good = all(cycleind1 ~= 0);
            for j = 1:4,
                cycleevent(cycleind1(j,good) + a-1) = j;
            end;
        end;
    end;
    
    freq = deriv(stimt,phase);

    if (istrans)
        phase = phase';
        freq = freq';
        amplitude = amplitude';
        cycleevent = cycleevent';
    end
    
    varargout = {phase,freq, amplitude, cycleevent};
  
  otherwise,
    error('Unrecognized phase estimation method');
end;

% -------------------------------------------------------------------------
% peaks and zeros function
function [phase,cycleind] = peaksandzeros(stimt,stim, smoothparam)

%find the peaks, anticipating that some peaks may be flattened by clipping
nflat = 1 / (stimt(2) - stimt(1));      % 1 sec
[hival,hi] = findpeaks2(stim,'max', 'numneighbors',3, 'maxflat',nflat);
[loval,lo] = findpeaks2(stim,'min', 'numneighbors',3, 'maxflat',nflat);

npeak = min(length(hival),length(loval));
med = mean((hival(1:npeak) + loval(1:npeak))/2);

%also find the zeros
[zc,sgn] = findzeros(stim-med);

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
    phase = NaN(size(stimt));
    cycleind = zeros(size(stimt));
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
            if (abs(stim(pts1(i))) > abs(stim(pts1(i+1))))
                good1(i+1) = false;
            else
                good1(i) = false;
            end;
        elseif (((tp1(i) == 1) || (tp1(i) == 3)) && ...
            ((tp1(i+1) == 1) || (tp1(i+1) == 3))),
            %keep the closer to zero of two repeated zeros
            if (abs(stim(pts1(i))) < abs(stim(pts1(i+1))))
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
cyclet(ind) = stimt(pts(good));
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
    good2 = (stimt >= min(cyclet(:,1))) & (stimt <= max(cyclet(:,end)));
    
    %interpolate the phase values
    phase = NaN(size(stimt));
    if (sum(good1(:)) > 3)
        h = nanmean(diff(cyclet(good1)));
        sp = csaps(cyclet(good1),ph1(good1), 1/(1+h^3/smoothparam));
        phase(good2) = fnval(sp, stimt(good2));
    end
else
    phase = NaN(size(stimt));
    cycleind = zeros(size(stimt));
end;





