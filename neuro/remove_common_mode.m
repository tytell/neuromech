function [rem,noise,noisemag] = remove_common_mode(sig, fs,goodband, varargin)

opt.order = 9;
opt.clip = 10;
opt = parsevarargin(opt, varargin, 2);

[z,p,k] = butter(opt.order,goodband/(fs/2), 'stop');
[sos,g] = zp2sos(z,p,k);

nchan = size(sig,2);
noise = filtfilt(sos,g, sig);

isclip = abs(sig) >= opt.clip;
noise(isclip) = NaN;

noisemag = sqrt(nanmean(noise.^2));

noise = noise ./ repmat(noisemag, [size(sig,1) 1]);

noisesign = ones(1,nchan);
for i = 2:nchan
    noisesign(i) = sign(nansum(noise(:,1).*noise(:,i)));
end
noisemag = noisemag.*noisesign;
noise = bsxfun(@times,noise,noisesign);

noise = nanmean(noise,2);

rem = sig - noise*noisemag;
