function varargout = wavecoherence(scales,coefs1,coefs2,wa,wb, varargin)
% function [coherence,phase] = wavecoherence(scales,coefs1,coefs2,wa,wb)
% Calculates coherence for two wavelet spectra and returns the squared
% coherence and the phase of the cross spectrum.  wa and wb are the
% smoothing parameters across scale and time, respectively

opt.mincoef = 0;
opt = parsevarargin(opt,varargin, 6);

W12 = coefs1 .* conj(coefs2);
W12 = avgwavespec(scales,bsxfun(@rdivide, W12,scales), wa,wb);

%we don't use coefs1 .* conj(coefs1) because of a weird issue with NaNs in
%complex numbers: NaN+NaNi .* conj(NaN+NaNi) = NaN + NaNi, which is still
%complex, but abs(NaN+NaNi)^2 = NaN, which is real
W11 = abs(coefs1).^2;
W11 = avgwavespec(scales,bsxfun(@rdivide, W11,scales), wa,wb);

W22 = abs(coefs2).^2;
W22 = avgwavespec(scales,bsxfun(@rdivide, W22,scales), wa,wb);

coherence = abs(W12).^2 ./ (W11 .* W22);

%we have problems when some coefficients are exactly zero.  But the
%coherence there is also zero, so just set them to zero
bad = (W11 <= opt.mincoef) | (W22 <= opt.mincoef);
coherence(bad) = 0;

varargout{1} = coherence;
if (nargout == 2),
    varargout{2} = angle(W12);
end;




