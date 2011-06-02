function varargout = waveletslide(s,dt,varargin)
% Continuous wavelet transform using a sliding FFT.  Much faster than how 
% the Wavelet toolbox implements the CWT, and somewhat faster than Torrence
% and Compo.  
% 
% Most of the code here taken from Torrence and Compo:
%    ``Wavelet software was provided by C. Torrence and G. Compo,
%      and is available at URL: http://paos.colorado.edu/research/wavelets/''.
%
% Reference: Torrence, C. and G. P. Compo, 1998: A Practical Guide to
%            Wavelet Analysis. <I>Bull. Amer. Meteor. Soc.</I>, 79, 61-78.

if (nargin == 2),
    error('Too few arguments');
end;

%defaults for the options
wavename = '';
waveparam = [];

ispad = true;
padlength = [];
iswindow = true;
isautoscale = true;
dj = [];
s0 = [];
j1 = [];
scales = [];
isdemean = true;
isdiagnostics = true;
filtersize = 6;

if (isnumeric(varargin{1})),
    scales = varargin{1};
    isautoscale = false;
    i = 2;
else
    i = 1;
end;

%options
while (i <= length(varargin)),
    switch lower(varargin{i}),
        case {'morlet','paul','dog'},
            wavename = varargin{i};
            if ((length(varargin) > i) && isnumeric(varargin{i+1}) && ...
                    (numel(varargin{i+1}) == 1)),
                waveparam = varargin{i+1};
                i = i+2;
            end;
        case 'pad',
            if (length(varargin) > i),
                if (islogical(varargin{i+1})),
                    ispad = varargin{i+1};
                    i = i+2;
                elseif (isnumeric(varargin{i+1})),
                    if ((numel(varargin{i+1}) ~= 1) || ...
                            (floor(varargin{i+1}) ~= varargin{i+1})),
                        error('Incorrect specification for pad/window length');
                    end;
                    ispad = true;
                    padlength = varargin{i+1};
                    i = i+2;
                elseif (ischar(varargin{i+1})),
                    switch lower(varargin{i+1}),
                        case {'on','true'},
                            ispad = true;
                            i = i+2;
                        case {'off','false'},
                            ispad = false;
                            i = i+2;
                        otherwise,
                            ispad = true;
                            i = i+1;
                    end;
                end;
            else
                ispad = true;
                i = i+1;
            end;

        case 'dt',
            dt = varargin{i+1};
            i = i+2;
        case {'dj','scalespacing'},
            dj = varargin{i+1};
            i = i+2;
        case 'nvoice',
            dj = 1/varargin{i+1};
            i = i+2;
        case {'s0','minscale'},
            s0 = varargin{i+1};
            i = i+2;
        case 'j1',
            j1 = varargin{i+1};
            i = i+2;
        case 'filtersize',
            filtersize = varargin{i+1};
            i = i+2;
        case 'numscales',
            j1 = varargin{i+1}-1;
            i = i+2;
        case {'autoscale','demean','window','diagnostics'},
            param = varargin{i};
            if (length(varargin) > i),
                if (islogical(varargin{i+1})),
                    val = varargin{i+1};
                    i = i+2;
                elseif (ischar(varargin{i+1})),
                    switch lower(varargin{i+1}),
                        case {'on','true'},
                            val = true;
                            i = i+2;
                        case {'off','false'},
                            val = false;
                            i = i+2;
                        otherwise,
                            val = true;
                            i = i+1;
                    end;
                end;
            else
                val = true;
                i = i+1;
            end;
            eval(['is' param '= val;']);

        otherwise,
            error('Unrecognized option %s\n', varargin{i});
    end;
end;

if (isempty(wavename)),
    error('Must specify a wavelet');
end;

N = length(s);

%remove the mean
if (isdemean),
    m = mean(s);
    s = s - m;
    dprintf(isdiagnostics,'Removing mean %g\n',m);    
end;

%pad the time series, if necessary
if (~isempty(padlength)),
    s(end+1:padlength) = 0;
elseif (ispad),
    padlength = 2^nextpow2(N);   % power of 2 closest to N
    s(end+1:padlength) = 0;
    dprintf(isdiagnostics,'Auto-padding to length %d\n',padlength);    
end;

%set up the scales if necessary
if (isautoscale),
    if (~isempty(scales)),
        warning('cwtfft:autoscaleconflict','Autoscaling option will overwrite scales parameter');
    end;
    if (isempty(s0)),
        s0 = 2*dt;
    end;
    if (isempty(dj)),
        dj = 0.25;
    end;
    if (isempty(j1)),
        j1 = fix((log(N*dt/s0)/log(2))/dj);
    end;
    
    %note, Matlab uses nondimensional scales (~samples) while T&C use
    %dimensional scales (~time).  We'll use dimensional scales, since we're
    %mostly following T&C
    scales = s0*2.^((0:j1)*dj);
    if (isdiagnostics),
        if (j1 < 30),
            fprintf('Using %d scales: ',j1);
            fprintf('%g ',scales);
            fprintf('\n');
        else
            fprintf('Using %d scales from %g to %g\n', j1,s0,scales(end));
        end;
    end;
    nscales = length(scales);
else
    nscales = length(scales);
    if (nscales == 0),
        error('If autoscaling is off, the scales must be passed as a parameter');
    end;
end;

%reset N
N = length(s);

coefs = complex(zeros(nscales,N));
blocksize2 = zeros(nscales,1);

islargescale = false;
for i = 1:length(scales),
    a = scales(i);
    
    nfilt = round(2 * a * sqrt(2*filtersize) / dt);
    
    if (~iswindow || (nfilt > N/2)),
        if (~islargescale),
            %wavenumber parameter
            k = (1:fix(N/2)) * 2*pi/(N*dt);
            k = [0 k -k(fix((N-1)/2):-1:1)];

            S = fft(s);
            islargescale = true;
        end;
        
        wavefft = getwaveletfft(wavename,k,a,waveparam);
        W = ifft(S.*wavefft);
        coefs(i,:) = W;
    else
        blocksize = estblocksize(N,nfilt);   
        if (blocksize > N+nfilt-1),
            blocksize = N+nfilt-1;
        end;
        blocksize2(i) = blocksize;
        L = blocksize - nfilt + 1;

        %wavenumber parameter
        k = (1:fix(blocksize/2)) * 2*pi/(blocksize*dt);
        k = [0 k -k(fix((blocksize-1)/2):-1:1)];

        wavefft = getwaveletfft(wavename,k,a,waveparam);

        %code below more or less duplicated from fftfilt.m
        istart = 1;
        while (istart <= N),
            iend = min(istart+L-1, N);
            if (iend - istart == 0),
                S = s(istart(ones(blocksize,1)));
            else
                S = fft(s(istart:iend),blocksize);
            end;

            W = ifft(S.*wavefft);
            wend = min(N,istart+blocksize-1);
            n1 = wend-istart+1;

            coefs(i,istart:wend) = coefs(i,istart:wend) + W(1:n1);

            istart = istart + L;
        end;
    end;
end;

switch wavename,
    case 'morlet',
        k0 = waveparam;
        fourier_factor = (4*pi)/(k0 + sqrt(2 + k0^2)); % Scale-->Fourier [Sec.3h]
        coi = fourier_factor/sqrt(2);                  % Cone-of-influence [Sec.3g]
end;

per = scales*fourier_factor;

%handle the outputs
switch nargout,
    case 1,
        varargout = {coefs};
        
    case 2,
        varargout = {coefs,per};

    case 3,
        varargout = {coefs,per,scales};
        
    case 4,
        varargout = {coefs,per,scales,blocksize2};
        
    otherwise,
        error('Incorrect number of outputs');
end;


function wavefft = getwaveletfft(wavename,k,scale,param)

n = length(k);
wavefft = zeros(1,n);

switch wavename, 
    case 'morlet',
        k0 = param;
        
        heaviside = k > 0;
        expnt = -(scale.*k(heaviside) - k0).^2 / 2;
        norm = sqrt(scale*k(2))*(pi^-0.25)*sqrt(n);     %normalization factor
        wavefft(heaviside) = norm*exp(expnt);
        
    case 'paul',
        m = param;
        
        heaviside = k > 0;
        expnt = -scale * k(heaviside);
        norm = sqrt(scale*k(2))*(2^m/sqrt(m*prod(2:(2*m-1))))*sqrt(n);
        wavefft(heaviside) = norm*exp(expnt);
        
    case 'dog',
        m = param;
        
        expnt = -(scale*k).^2 / 2;
        norm = sqrt(scale*k(2) / gamma(m+0.5))*sqrt(n);
        wavefft = -norm*(i^m)*((scale*k).^m).*exp(expnt);
end;


function blocksize = estblocksize(N,nfilt)
%from fftfilt.m

fftflops = [ 18 59 138 303 660 1441 3150 6875 14952 32373 69762 ...
    149647 319644 680105 1441974 3047619 6422736 13500637 28311786 59244791];

n = 2.^(1:20);
validset = n > (nfilt-1);   % must have nfft > (nb-1)
n = n(validset);
fftflops = fftflops(validset);

% minimize (number of blocks) * (number of flops per fft)
L = n - (nfilt - 1);
[dum,ind] = min( ceil(N./L) .* fftflops ); %#ok
blocksize = n(ind);


function dprintf(isdiagnostics,tplt,varargin)

if (isdiagnostics),
    fprintf(tplt,varargin{:});
end;


