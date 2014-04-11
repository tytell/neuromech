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
            [phase1,cycleind1] = getphase(stimt(istreat),stim1-med, 'smoothparam',opt.smoothparam);

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


