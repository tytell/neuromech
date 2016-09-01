function F = frequencypower(t,sig, freqs)

F = zeros(size(freqs));
nfreq = length(freqs);
N = length(t);

fsig = complex(zeros(N,nfreq));
for i = 1:length(freqs)
    %complex exponential at frequency freqs(i)
    fsig(:,i) = exp(-2*pi*1i*freqs(i));
end

F = sum(repmat(sig,[1 nfreq]) .* fsig);



    