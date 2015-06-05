function siglo = get_low_baseline(t,sig, locut)
% GET_LOW_BASELINE
% Gets a low frequency baseline for a signal by taking an average in
% blocks, then constructing a smooth signal

dt = t(2) - t(1);
lenwind = round(1/locut/dt);

if (size(sig,1) < size(sig,2))
    warning('Signal should be arranged in columns. Results may be weird')
end

if lenwind > size(sig,1)
    siglo = repmat(nanmean(sig),[size(sig,1) 1]);
    return;
end

%pad signal to get the length evenly divisible by lenwind
nextra = ceil(size(sig,1)/lenwind)*lenwind - size(sig,1);

npre = floor(nextra/2);
npost = ceil(nextra/2);

preval = nanmean(sig(1:lenwind,:));
postval = nanmean(sig(end-lenwind+1:end,:));

sigpad = cat(1, repmat(preval,[npre 1]), ...
    sig, ...
    repmat(postval,[npost 1]));

nwind = size(sigpad,1) / lenwind;
sz = size(sigpad);
sz = [lenwind nwind sz(2:end)];

sigpad = reshape(sigpad,sz);

sigmn = nanmean(sigpad,1);
sigmn = squeeze(sigmn);
if ((size(sigmn,1) == 1) && (size(sigmn,2) > 1))
    sigmn = sigmn';
end

sigmn = sigmn([1 1:end end],:);
ctrind = cat(1,0,(1:nwind)'*lenwind - lenwind/2 - npre, size(sig,1)+1);
ind = (1:size(sig))';

siglo = interp1(ctrind,sigmn, ind, 'spline');




