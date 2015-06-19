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

if opt.showclosefreqs
    for a = 1:length(freq)
        freq1 = freq{a};
        fprintf('** ');
        fprintf('%f ', freq1);
        fprintf('\n');
        
        m = abs(repmat(freq1(:),[1 length(freq1)]) - repmat(freq1(:)',[length(freq1) 1])) / basefreq * Nsec;
        p = repmat(freq1(:),[1 length(freq1)]) + repmat(freq1(:)',[length(freq1) 1]) / basefreq * Nsec;
        
        %set the diagonal to Inf
        isdiag = eye(length(freq1)) == 1;
        m(isdiag) = Inf;

        m = mod(m + 0.5, 1) - 0.5;
        p = mod(p + 0.5, 1) - 0.5;
        
        [~,mord] = sort(abs(m(:)));
        [closei,closej] = ind2sub(size(m),mord);
        
        fprintf('10 closest differences\n');
        for i = 1:10
            fprintf('  %f, %f: %f\n', freq1(closei(i)), freq1(closej(i)), m(mord(i)));
        end
        
        [~,pord] = sort(abs(p(:)));
        [closei,closej] = ind2sub(size(p),pord);
        
        fprintf('10 closest sums\n');
        for i = 1:10
            fprintf('  %f, %f: %f\n', freq1(closei(i)), freq1(closej(i)), p(pord(i)));
        end
    end
end
