function varargout = findbursts(t,varargin)
%Finds bursts based a series of spike times (and optionally spike heights,
%assuming that noise is smaller than spikes).  Uses a variety of
%techniques:
%   'runmean' - Running mean spike time technique.  Looks for regions in
%   which the running mean spike time changes slowly.  Takes one parameter,
%   the running mean duration, which should be approximately the cycle
%   period
%   'surprise' - Poisson surprise technique
%   'timing' - Looks for gaps in spike times, collecting bursts based on
%   minimum interburst intervals, etc.
%   'integral' - Looks for when the integrated neural signal passes above
%   or below a threshold
%   'peaks' - Looks for peaks in the spike rate

opt.isWeighted = false;

%if they bother to pass in the spike heights, assume they want means
%weighted by spike height
if ((nargin >= 2) && (isnumeric(varargin{1}) || iscell(varargin{1}))),
    spikeh = varargin{1};
    varargin = varargin(2:end);
    opt.isWeighted = true;
else
    spikeh = [];
end;

opt.stepFraction = 0.1;
opt.dutycyle = 0.3;
opt.cutoff = 0.2;
opt.minspikes = 0;
opt.cycledur = [];
opt.rundur = 1;
opt.correctmissedbursts = false;
opt.threshold = [];
opt.mergegap = 0;
opt.binsize = 0.05;
opt.smooth = 0.05;
opt.maxflat = 3;
opt.robustonoff = [];
opt.isoutstruct = true;
opt.isoutcell = false;
opt.mergemultiples = true;
opt.interburstdur = [];
opt.interburstaction = 'merge';

opt.method = 'runmean';

%check for options
p = 1;
while (p <= length(varargin)),
    switch lower(varargin{p}),
      case 'runmean',
        opt.method = varargin{p};
        opt.rundur = varargin{p+1};
        p = p+2;
        
      case {'poisson','surprise','timing','integral','peaks'},
        opt.method = varargin{p};
        p = p+1;
        
      case {'interburstdur','burstdur','interburstaction'},
        opt.(lower(varargin{p})) = varargin{p+1};
        p = p+2;
        
      case {'burstfreq','dutycycle','cutoff','minspikes','cycledur'},
        opt.(lower(varargin{p})) = varargin{p+1};
        p = p+2;
        
      case 'threshold',
        opt.threshold = varargin{p+1};
        p = p+2;

      case 'mergebursts',
        opt.mergegap = varargin{p+1};
        p = p+2;
        
      case 'correctmissedbursts',
        opt.correctmissedbursts = true;
        p = p+1;
        
        case {'outstruct','outcell'},
            opt.(['is' lower(varargin{p})]) = true;
            p = p+1;
            
      case 'weight',
        switch lower(varargin{p+1}),
          case 'on',
            opt.isWeighted = true;
          case 'off',
            opt.isWeighted = false;
          otherwise,
            error('Unrecognized option for weight %s.',varargin{p+1});
        end;
        p = p+2;
    
      case {'smooth','binsize','maxflat','robustonoff'},
        opt.(lower(varargin{p})) = varargin{p+1};
        p = p+2;
        
      case {'mergemultiples','israte'},
        opt.(lower(varargin{p})) = true;
        p = p+1;
        
      otherwise,
        error('Unrecognized option %s.',varargin{p});
    end;
end;

if (nargout ~= 2)
    opt.isoutstruct = false;
end;

if (length(t) < 10),
    warning('Too few spikes for burst detection.');
    varargout = cell(nargout,1);
    return;
end;

%check for irregular spacing in the time, which indicates that we have a
%spike train
dt = t(2)-t(1);
if (any(abs((diff(t)-dt)/dt) > 0.01)),
    opt.isSpikes = true;
else
    opt.isSpikes = false;
end;

switch opt.method,
  case 'runmean',
    [burstctr burstdev burstskew spikectr tctr] = findbursts_runmean(t,spikeh,opt);
        
    if (nargout == 1),
        varargout = {burstctr};
    elseif (nargout == 3),
        varargout = {burstctr,burstdev,burstskew};
    elseif (nargout == 5),
        varargout = {burstctr,burstdev,burstskew,spikectr,tctr};
    end;
        
  case 'surprise',
    [burstctr, burstind, burstSurprise, surprise] = findbursts_surprise(t,spikeh, opt);
    varargout = {burstctr, burstind, burstSurprise, surprise};
    
  case {'timing','integral'},
    [burstctr, burstind] = findbursts_timing(t,spikeh, opt);
    
    varargout = {burstctr,burstind};

  case 'peaks',
      if (numel(t) == size(t,1)),
          t = t';
          spikeh = spikeh';
          nchan = 1;
      else
          nchan = size(t,1);
      end;
      
      for ch = 1:nchan,
          good = isfinite(t(ch,:));
          if (~isempty(spikeh))
              good = good & isfinite(spikeh(ch,:));
              spikeh1 = spikeh(ch,good);
          else
              spikeh1 = [];
          end;
          t1 = t(ch,good);
          
          [burstctr{ch}, burstind{ch}, peakind{ch},peakrate{ch},tbin{ch},...
              spikerate{ch},ismerged{ch}] = ...
              findbursts_peaks(t1,spikeh1, opt);
      end;
      
      if (opt.isoutstruct),
          spike.tbin = tbin;
          spike.rate = spikerate;
          
          burst.ctr = burstctr;
          burst.spikeratepk = peakrate;
          
          for ch = 1:nchan,
              good = isfinite(t(ch,:));
              spike.t{ch} = t(ch,good);
              if (~isempty(spikeh)),
                  spike.h{ch} = spikeh(ch,good);
              end;

              burst.on{ch} = t(burstind{ch}(1,:));
              burst.off{ch} = t(burstind{ch}(2,:));

              n = size(burstind{ch},2);
              
              %number of spikes in the burst divided by the duration gives us the mean spike rate
              burst.spikeratemn{ch} = (diff(burstind{ch})+1) ./ ...
                  (burst.off{ch} - burst.on{ch});
              
              spike.burstind{ch} = zeros(1,size(spike.t,2));
              for b = 1:n,
                  k = burstind{ch}(1,b):burstind{ch}(2,b);
                  
                  %set the burst index for each spike
                  spike.burstind{ch}(k) = b;
              end;
          end;
          
          if (~opt.isoutcell),
              burst = structfun(@(x) (catuneven(1,x{:})), burst,'UniformOutput',false);
              spike = structfun(@(x) (catuneven(1,x{:})), spike,'UniformOutput',false);              
          elseif (nchan == 1),
              burst = structfun(@(x) (x{1}), burst);
              spike = structfun(@(x) (x{1}), spike);
          end;
          varargout = {burst,spike};
      else
          if (nchan == 1),
              burstctr = burstctr{1};
              burstind = burstind{1};
              peakind = peakind{1};
              peakrate = peakrate{1};
              tbin = tbin{1};
              spikerate = spikerate{1};
              ismerged = ismerged{1};
          end;
          
          varargout = {burstctr,burstind,peakind,peakrate,tbin,spikerate,ismerged};
          varargout = varargout(1:nargout);
      end;
      
  case 'poisson',
    error('Not implemented yet.');
end;

            
