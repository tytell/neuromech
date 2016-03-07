function freq = get_pert_frequencies(basefreq,freqrng,ncycles,Nsec, varargin)
%GET_PERT_FREQUENCIES    Get perturbation frequencies for Ankarali code
%   freq = get_pert_frequencies(basefreq,freqrng,ncycles,Nsec)
%
%       basefreq = Base oscillation frequency in Hz
%       freqrng = Range of perturbation frequencies to identify.
%                 e.g. [1 20] means to identify from 1Hz to 20Hz.
%                 If the first element is zero, it will identify down to
%                 the lowest multiple of the frequency resolution
%       ncycles = Number of base oscillation cycles that will be run
%       Nsec = Number of Poincare sections to use in the cycle
%           If Nsec is NaN, then it uses prime multiples of the minimum
%           frequency
%
%   returns freq: a cell array of different combinations of frequenciues to
%    try.
%
% EXAMPLE
%   >> freq = get_pert_frequencies(3,[1 15],5,6)
%
% freq =
% 
%   [  1.8000    9.6000   13.8000  ]
%   [  4.8000    5.4000    8.4000   10.8000   ]
%   [ 10.2000   12.6000  ]
%   [ 13.2000   14.4000  ]
%   [  1.2000    2.4000    4.2000   ]
%   [  7.8000   11.4000   ]
%   [  6.6000    7.2000   ] 
%   [  3.6000  ]
%
% So we should run 8 trials of 5 cycles of oscillation at 3Hz with
% different combination of perturbation frequencies on top.  The first
% trial should have perturbations at 1.8, 9.6, and 13.8 Hz.

opt.randomize = true;
opt.keepfreq = 1;
opt.tol = 0.01;
opt.showclosefreqs = true;
opt.nclosefreqs = 10;
opt.makesignal = true;
opt.sampfreq = 10000;
opt.pertamp = 0.05;

opt = parsevarargin(opt,varargin, 5);

if length(freqrng) > 2
    freq1 = freqrng / basefreq;
    Nsec = 2;
elseif isnan(Nsec)
    fres1 = 1/ncycles;
    
    Nsec = 2;
    
    lofreq = freqrng(1)/2;
    assert(lofreq ~= 0);
    
    mult = ceil(freqrng(2) / lofreq);
    primemult = primes(mult);
    
    freq1 = (lofreq * primemult) / basefreq;
else    
    %get frequency resolution
    fres1 = 1/ncycles;
    %convert max freq
    maxfreq1 = freqrng(2) / basefreq;
    
    %make sure the minimum frequency is not zero and is a multiple of fres1
    minfreq1 = ceil(freqrng(1)/basefreq/fres1)*fres1;
    if minfreq1 < fres1
        minfreq1 = fres1;
    end
    %here's the range of nondimensional frequencies that we'll identify
    freq1 = minfreq1:fres1:maxfreq1;
end

if opt.keepfreq ~= 1
    nkeep = floor(opt.keepfreq * length(freq1));
    indkeep = round(linspace(1,length(freq1),nkeep));

    good = [true diff(indkeep) ~= 0];
    indkeep = indkeep(good);

    freq1 = freq1(indkeep);
end

%we can't identify multiples of 0.5, so get rid of them
isintfreq1 = divisible(freq1, 0.5);
freq1 = freq1(~isintfreq1);

%randomize the order
if opt.randomize
    ord = randperm(length(freq1));
    freq1 = freq1(ord);
end

%now look for groups of frequencies where the differences and sums of the
%frequencies, multiplied by Nsec, are not integers
freq = {};
a = 1;
while ~isempty(freq1)
    good = true(size(freq1));
    for i = 1:length(freq1)
        if good(i)
            %eliminate frequencies that don't meet the criteria relative to
            %the current one
            m = abs(freq1(i) - freq1(good))*Nsec;
            p = (freq1(i) + freq1(good))*Nsec;

            ism = false(size(freq1));
            isp = false(size(freq1));
            ism(good) = abs(m - round(m)) < opt.tol;
            isp(good) = abs(p - round(p)) < opt.tol;

            good(ism | isp) = false;
            good(i) = true;
        end
    end

    freq{a} = freq1(good);
    a = a+1;

    freq1 = freq1(~good);
end
    
for a = 1:length(freq)
    freq{a} = freq{a} * basefreq;
end

for a = 1:length(freq)
    freq1 = freq{a};
    fprintf('** ');
    fprintf('%f ', freq1);
    fprintf('\n');

    if length(freq1) == 1
        continue;
    end

    if opt.showclosefreqs
        freq1 = freq1 / basefreq;
        
        m = abs(repmat(freq1(:),[1 length(freq1)]) - repmat(freq1(:)',[length(freq1) 1])) * Nsec;
        p = (repmat(freq1(:),[1 length(freq1)]) + repmat(freq1(:)',[length(freq1) 1])) * Nsec;
        
        %get rid of the lower triangular part of the matrix
        for i = 1:length(freq1)
            for j = 1:length(freq1)
                if i <= j
                    m(i,j) = Inf;
                    p(i,j) = Inf;
                end
            end
        end
        
        m = abs(m - round(m));
        p = abs(p - round(p));
        
        [~,mord] = sort(m(:));
        [closei,closej] = ind2sub(size(m),mord);
        
        nclose = min(length(freq1),opt.nclosefreqs);
        fprintf('%d closest differences\n', nclose)
        for i = 1:nclose
            fprintf('  %f, %f: %f\n', freq1(closei(i))*basefreq, freq1(closej(i))*basefreq, m(mord(i)));
        end
        
        [~,pord] = sort(p(:));
        [closei,closej] = ind2sub(size(p),pord);
        
        fprintf('%d closest sums\n',nclose);
        for i = 1:nclose
            fprintf('  %f, %f: %f\n', freq1(closei(i))*basefreq, freq1(closej(i))*basefreq, p(pord(i)));
        end
    end
end

if opt.makesignal
    dt = 1/opt.sampfreq;
    
    t = (0:dt:ncycles/basefreq)';

    N = length(t);
    
    for a = 1:length(freq)
        sig1 = sin(2*pi*basefreq*t);
        
        freq1 = freq{a};
        ph = rand(length(freq1),1) * 2*pi;
        
        for i = 1:length(freq1)
            pert1 = opt.pertamp * sin(2*pi*freq1(i)*t + ph(i));
            sig1 = sig1 + pert1;
        end
        
        y = fft(sig1);
        Py = y .* conj(y) / N;

        f = 1/dt * (0:ceil(N/2))'/N;
        Py = Py(1:length(f),:);
        
        subplot(length(freq),1,a);
        plot(f,10*log10(Py));
        
        hb = vertplot(basefreq, 'r-');
        hp = vertplot(freq1,'k--');
        
        [~,k] = min(abs(f-basefreq));
        
        axis([0 1.2*freqrng(2) -20 10*log10(Py(k))]);
        xlabel('Frequency');
        ylabel('Power (dB)');
        
        if (a == length(freq))
            legend([hb hp(1)], 'Base','Perturbations');
        end
    end
end

    
    