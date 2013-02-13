function [coeffs,sampfreq,lofreqs,scales] = wavelofreqs(S,sampfreq,freqs,wavename)

decompwavelet = 'db3';
freqmultiple = 8;

fac = sampfreq / (freqmultiple * max(freqs));
Ndecomp = floor(log2(fac));

[C,L] = wavedec(S,Ndecomp,decompwavelet);
A = appcoef(C,L,decompwavelet);

sampfreq = sampfreq / (2^Ndecomp);

[scales,lofreqs] = frq2scal(freqs, wavename, 1/sampfreq);
[lofreqs,ord] = unique(lofreqs);
scales = scales(ord);

coeffs = cwt(A,scales, wavename);










