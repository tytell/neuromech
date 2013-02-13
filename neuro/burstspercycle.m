function varargout = burstspercycle(burstt, stimt,stim, varargin) 
% function [nburst,entrainlength,cyclestart,phase] = burstspercycle(burstt, stimt,stim, varargin)
%  or      [nburst,nburstadj,entrainlength,cyclestart,phase] = burstspercycle(...,'burstint',b)

opt.phase = [];
opt.phaseestimmode = 'phaser';
opt.smoothphase = true;
opt.appxcycleper = 1;
opt.burstint = [];

opt = parsevarargin(opt, varargin, 3, ...
                    'multival',{'phaseestimmode',{'phaser','peaks'}});

stimt = makerow(stimt);
stim = makerow(stim);

if (~isempty(opt.phase)),
    phase = opt.phase;
else
    switch lower(opt.phaseestimmode),
      case 'phaser',
        phr = newPhaser(stim);
        phase = feval(phr.phi, phr, stim);
        phase = phase / (2*pi);
        
      case 'peaks',
        error('Not implemented yet');
    end;
end;
phase = makerow(phase);

dt = stimt(2) - stimt(1);
sampfreq = 1/dt;

phase = phase - min(phase);
if (opt.smoothphase),
    n = round(0.1*opt.appxcycleper * sampfreq);
    phase = runavg(phase,n,2);
end;

cyclenum = floor(phase) + 1;
ncycle = max(cyclenum);

good = isfinite(cyclenum);
cyclestart = accumarray(cyclenum(good)',stimt(good)', [], @min);

cycleend = accumarray(cyclenum(good)',stimt(good)', [], @max);
goodcycle = accumarray(cyclenum(good)',stimt(good)', [], @(x) (range(x) == length(x)*dt));

burstphase = NaN(size(burstt));
good1 = isfinite(phase);
good2 = isfinite(burstt);

burstphase(good2) = interp1(stimt(good1),phase(good1), burstt(good2));

nburst = zeros(size(burstt,1),length(cyclestart));
if (~isempty(opt.burstint))
    nburstadj = zeros(size(burstt,1),length(cyclestart));
end;
for ch = 1:size(burstt,1),
    good = isfinite(burstphase(ch,:));
    nburst1 = accumarray(floor(burstphase(ch,good))'+1, 1, [ncycle 1]);
    nburst(ch,:) = nburst1;
    
    if (~isempty(opt.burstint)),
        nburstadj(ch,:) = accumarray(floor(burstphase(ch,good))'+1, opt.burstint(ch,good)', ...
                                     [ncycle 1], @adjnburst);
    end;
end;

entrainlength = zeros(size(nburst));
run1 = nburst == 1;
for ch = 1:size(nburst,1),
    a = 1;
    while (a <= size(nburst,2)),
        if (run1(ch,a)),
            b = a+1;
            while ((b <= size(nburst,2)) && run1(ch,b)),
                b = b+1;
            end;
            entrainlength(ch,a:b-1) = b-a;
            
            a = b;
        else
            a = a+1;
        end;
    end;   
end;

if (size(stimt,1) == 1),
    cyclestart = makerow(cyclestart);
end;

if (isempty(opt.burstint)),
    varargout = {nburst, entrainlength, cyclestart, phase};
else
    varargout = {nburst, nburstadj, entrainlength, cyclestart, phase};
end;

%------------------------------------------
function f = adjnburst(burstint)

if (isempty(burstint)),
    f = 0;
elseif (length(burstint) == 1)
    f = Inf;
else
    burstint = sort(burstint,'descend');
    f = burstint(1)/burstint(2);
end;

    
    
    



