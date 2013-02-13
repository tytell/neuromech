function [per,perbychan,chanstd] = robustper(t,varargin)
% function [per, perbychan,chanstd] = robustper(t,options...)
% t is a vector or matrix of events (in columns).  robustper estimates the average period
% between events by calculating increasingly higher order differences (t(2)-t(1) is first
% order, t(3)-t(1) is second, etc) and searching for a value of the period such that most
% of the differences are integer multiples.
%
% Options:
% 'maxdiff' - maximum order difference to use.  Default is as many as are possible based
% on the length of the data set.
% 'guessfrac' - Range of fractions of the first order difference to guess as the true
% period estimate.  Default is 0.8 to 1.3 in steps of 0.05.
% 'channeldim' - Dimension that represents multiple channels, sampled simultaneously, that
% ought to have the same period.  Merges the channels to produce a joint estimate of
% average period.

opt.maxdiff = Inf;
opt.guessfrac = 0.8:0.05:1.3;
opt.channeldim = [];
opt.weightchannels = true;

i = 1;
while (i <= length(varargin)),
    switch lower(varargin{i}),
      case {'maxdiff','guessfrac','channeldim'},
        opt.(varargin{i}) = varargin{i+1};
        i = i+2;
        
      case {'weightchannels'},
        if ((i+1 <= length(varargin)) && islogical(varargin{i+1})),
            opt.(varargin{i}) = varargin{i+1};
            i = i+2;
        else
            opt.(varargin{i}) = true;
            i = i+1;
        end;
                                
      otherwise,
        error('Unrecognized option %s', varargin{i});
    end;
end;

%turn t into a column, if necessary
if (size(t,1) == 1),
    t = t(:);
    istranspose = true;
else
    istranspose = false;
end;

n = size(t,1);
sz = size(t);
if (~isempty(opt.channeldim)),
    pmt = [1 opt.channeldim 2:opt.channeldim-1 opt.channeldim+1:ndims(t)];
    nchan = size(t,opt.channeldim);
else
    pmt = [1 ndims(t)+1 2:ndims(t)];
    nchan = 1;
    sz(end+1) = 1;
end;

t = permute(t,pmt);
sz = sz(pmt);
t = reshape(t,[n nchan prod(sz(3:end))]);
N = size(t,3);

% maximum possible difference
if (~isfinite(opt.maxdiff)),
    opt.maxdiff = n-1;
end;

%total number of differences, first order to nth
np = sum(n-opt.maxdiff:n-1);
P = zeros(np,nchan,N);                  % unscaled period estimates

a = 0;
for off = 1:np,
    p = t(1+off:end,:,:) - t(1:end-off,:,:);

    if (off == 1),
        %use the first order differences to estimate the first guess for the period
        p0 = nanmedian(flatten(p,1:2));
    end;
    
    P(a+(1:size(p,1)),:,:) = p;
    
    a = a+size(p,1);
end;

%look for the appropriate integer multiple for the channels jointly
P = reshape(P,[np*nchan N]);
p0 = p0(ones(np*nchan,1),:);

%now, estimate the integer multiple for each difference by dividing by the guess for
%the period
%if all of the periods, divided by the integer multiple, work out to something close to
%the original guess, then we've got a good estimate
err = zeros(length(opt.guessfrac),N);
for i = 1:length(opt.guessfrac),
    n = round(P ./ (opt.guessfrac(i) * p0));
    n(n == 0) = 1;
    
    err1 = abs((P./n)-(opt.guessfrac(i)*p0));
    err1(isnan(err1)) = 0;
    err(i,:) = sum(err1);
end;

%minimum error
[q,ind] = min(err);

%and the fraction of the first order guess that corresponded to the minimum error
frac = opt.guessfrac(ind);
p = p0 .* frac(ones(np*nchan,1),:);

%estimate the true period based on those integer multiples
n = round(P ./ p);
n(n == 0) = NaN;

%reshape so that we can look at the channels separately again
p = P ./ n;

p = reshape(p,[np nchan N]);
perbychan = nanmedian(p);
chanstd = nanstd(p);

if (opt.weightchannels),
    weight = 1./(chanstd.^2);
    per = nansum(perbychan.*weight,2) ./ nansum(weight,2);
else
    per = nanmedian(p);
end;

per = reshape(per,[1 1 sz(3:end)]);
per = per(:,ones(1,nchan),:);
per = ipermute(per,pmt);

perbychan = reshape(perbychan,[1 nchan sz(3:end)]);
perbychan = ipermute(perbychan,pmt);

chanstd = reshape(chanstd,[1 nchan sz(3:end)]);
chanstd = ipermute(chanstd,pmt);

if (istranspose),
    per = per';
    chanstd = chanstd';
    perbychan = perbychan';
end;



