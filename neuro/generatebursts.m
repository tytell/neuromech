function [t,sig,phase] = generatebursts(sampfreq,len, varargin)
% Generates sample bursts based on a set of parameters.  Each parameter can
% be a one or two element matrix.  The first element specifies the mean
% value, and the optional second one specifies a standard deviation.
%
% Parameters:
%    bursttimefreq - Two column matrix.  First column is time, second is
%    burst frequency.  No standard deviations are allowed.
%    burstper - Burst period
%    burstdur - Burst duration, specified in the same units as burst period
%    dutycycle - Burst duration in fractions of the cycle
%    spikeheight - Spike height
%    spikeper - Period between spikes, on average

opt.burstfreq = 1;
opt.burstdur = [0.3 0];
opt.spikeheight = [1 0];
opt.spikeper = [0.005 0];
opt.bursttimefreq = [];
opt.chanphase = [0 0; 0.5 0];

%process the parameters
param = 1;
while (param <= length(varargin)),
    switch lower(varargin{param}),
      case {'burstfreq','burstdur','spikeheight', ...
            'spikeper','nchan'},
        nm = lower(varargin{param});
        opt.(nm) = varargin{param+1};
        param = param+2;
        
      case {'bursttimefreq'},
        nm = lower(varargin{param});
        opt.(nm) = varargin{param+1};
        param = param+2;
        
      case {'chanphase',},
        nm = lower(varargin{param});
        opt.(nm) = varargin{param+1};
        param = param+2;

      otherwise,
        error('Unrecognized option %s', varargin{param});
    end;
end;

if (~isempty(opt.chanphase)),
    nchan = size(opt.chanphase,1);
else
    nchan = 1;
end;

if (size(opt.burstdur,1) == 1),
    opt.burstdur = opt.burstdur(ones(nchan,1),:);
end;
if (size(opt.spikeheight,1) == 1),
    opt.spikeheight = opt.spikeheight(ones(nchan,1),:);
end;
if (size(opt.spikeper,1) == 1),
    opt.spikeper = opt.spikeper(ones(nchan,1),:);
end;
    
if (isempty(opt.bursttimefreq)),
    opt.bursttimefreq = [0 opt.burstfreq; len opt.burstfreq];
end;

t = (0:round(len*sampfreq))/sampfreq;
sig = zeros(nchan,length(t));

%generate a continuous burst frequency
burstfreq = spline(opt.bursttimefreq(:,1),opt.bursttimefreq(:,2), t);

%and integrate to produce a phase
phase = cumtrapz(t,burstfreq);
    
%total number of bursts
nburst = floor(phase(end))-1;

burstphase = zeros(nchan,nburst);
for ch = 1:nchan,
    burstphase(ch,:) = opt.chanphase(ch,1) + randn(1,nburst) * opt.chanphase(ch,2) + ...
        (1:nburst);
end;

burstind = zeros(nchan,nburst);

wrap = find(diff(mod(phase,1)) < -0.5);
wrap = [1 wrap];

for b = 1:nburst,
    for ch = 1:nchan,
        burstind(ch,b) = first(phase >= burstphase(ch,b), wrap(b));
    end;
end;

%generate the burst lengths in number of samples
burstlen = abs(opt.burstdur(:,ones(1,nburst)) + randn(nchan,nburst) .* opt.burstdur(:,2*ones(1,nburst)));

%run through each burst
for i = 1:nburst,
    for ch = 1:nchan,
        if (opt.spikeper(1) > 0),
            %approximate number of spikes
            nspike = floor(burstlen(i)/opt.spikeper(ch,1));
            
            spikegap = abs(opt.spikeper(ch,1) + randn(1,nspike) * opt.spikeper(ch,2));
            spiket = cumsum(spikegap);
            spiket = spiket(spiket <= burstlen(ch,i));
            spikeind = round(spiket * sampfreq);
            
            spikeval = abs(opt.spikeheight(ch,1) + randn(1,length(spikeind)) * ...
                           opt.spikeheight(ch,2));
        else
            spikeind = 0:burstlen(ch,i)-1;
            spikeval = ones(size(spikeind));
        end;

        spikeind = spikeind + burstind(ch,i);
        
        if (i == nburst),
            spikeind = spikeind(spikeind <= length(t));
        end;
        sig(ch,spikeind) = spikeval;
    end;
end;


    
        
        
    
    


