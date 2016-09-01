function F = fouriertransform(t,sig, freqs)
%FOURIERTRANSFORM   Directly evaluate Fourier transform at specific
%frequncies
%
%  Does not use FFT.  Evaluates the Fourier transform directly at a set
%  of specific frequencies.

nfreq = length(freqs);
N = length(t);

fsig = complex(zeros(N,nfreq));
for i = 1:length(freqs)
    %complex exponential at frequency freqs(i)
    fsig(:,i) = 2/N * exp(-2*pi*1i*t*freqs(i));
end

sz = size(sig);
F = zeros(nfreq,sz(2:end));
nreps = prod(sz(2:end));
for i = 1:nreps
    sig1 = sig(:,i) - nanmean(sig(:,i));
    F1 = sum(repmat(sig1,[1 nfreq]) .* fsig);
    F(:,i) = F1';
end




    