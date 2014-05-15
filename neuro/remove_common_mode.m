function [rem,common,commonmag] = remove_common_mode(sig, varargin)

opt.nfft = 512;
opt.window = [];
opt.noverlap = [];
opt.coherethresh = 0.1;
opt = parsevarargin(opt, varargin, 2);

nchan = size(sig,2);
coh = zeros(opt.nfft,nchan*(nchan-1)/2);
k = 1;
for i = 1:nchan
    for j = i+1:nchan
        [coh(:,k),f1] = mscohere(sig(:,i),sig(:,j),hanning(opt.nfft),opt.noverlap,opt.nfft,'twosided');
        k = k+1;
    end
end

cohall = median(coh,2);
f2 = linspace(0,2*pi,size(sig,1)+1)';
f2 = f2(1:end-1);
cohall = interp1(f1,cohall, f2, 'linear',0);

F = fft(sig);
common = zeros(size(sig));
iscommon = cohall > opt.coherethresh;
commonsign = ones(1,nchan);
for i = 1:nchan
    common(:,i) = real(ifft(F(:,i).*double(iscommon)));
    if (i > 1)
        commonsign(i) = sign(sum(common(:,1).*common(:,i)));
    end
end
rem = sig - common;
