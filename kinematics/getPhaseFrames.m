function [phase,phaseframes,nextbeatfr] = getPhaseFrames(t,indpeak,tx, nphase)
% function [phase,phaseframes,nextbeatfr] = ...
%                        getPhaseFrames(t,indpeak,tx, nphase)

% only operate on peaks at the tail
indpeak = indpeak(end,:);
npeak = length(indpeak);

% check to see whether we start with the tail on the left or the right
% we arbitrarily define phase = 0 to be on the left side (x coordinate
% high)
if (tx(indpeak(1)) < tx(indpeak(2))),
    % we're starting on the right side, so initial phase = pi
    startphase = pi;
else
    startphase = 0;
end;

sp = csapi(indpeak(1:end-2),indpeak(3:end));
nextbeatfr = repmat(NaN,size(t));
nextbeatfr(indpeak(1):indpeak(end-2)) = ...
    fnval(sp,indpeak(1):indpeak(end-2));

phase = repmat(NaN,size(t));

% we use the times defined by indpeak as instances of known phase (starting
% at either 0 or pi, depending on which side the tail began on) and then
% increasing in increment of pi.  Then we can use spline to generate a
% smooth phase function
phase0 = startphase + (0:npeak-1)*pi;
phase(indpeak(1):indpeak(end)) = ...
    mod(spline(indpeak,phase0, indpeak(1):indpeak(end)),2*pi);

% now we need to invert the phase function - determine which frames
% correspond to a particular phase
% we use spline again

% first define the phase values we want to know about
pfinv = repmat(linspace(0,2*pi,nphase+1)',[1 ceil(npeak/2)+1]);
pfinv = pfinv + repmat(2*pi*(0:ceil(npeak/2)),[nphase+1 1]);
pfinv = pfinv(1:end-1,:);

% then use spline to estimate the frames corresponding to those phases
phaseframes = repmat(NaN,size(pfinv));
phaseframes = spline(phase0,indpeak, pfinv);
phaseframes = round(phaseframes);

% remove any extrapolations
phaseframes(pfinv > phase0(end)) = NaN;
phaseframes(pfinv < startphase) = NaN;

% eliminate any all NaN columns
k = find(any(isfinite(phaseframes)));
phaseframes = phaseframes(:,k);

