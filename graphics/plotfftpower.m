function varargout = plotfftpower(t,sig, varargin)
%PLOTFFTPOWER  Plots an FFT power spectrum
%
% [f,Py] = plotfftpower(t,sig)
%  or      plotfftpower(sig,sampfreq)
%  or      plotfftpower(sig)        (assumes sampling freq = 1Hz)
%

opt.scale = 'linear';       % or 'dB'

if ((nargin >= 2) && isnumeric(t) && isnumeric(sig)),
    if (length(sig) == 1),
        dt = 1/sig;
        sig = t;
        t = [];
    else
        dt = t(2)-t(1);
    end;
elseif (nargin == 1) || ((nargin >= 1) && ischar(sig)),
    varargin = [sig varargin];
    sig = t;
    dt = 1;
    t = [];
end;

if ~isempty(varargin) && matchlinespec(varargin{1})
    ls = varargin{1};
    varargin = varargin(2:end);
else
    ls = '-';
end

opt = parsevarargin(opt,varargin, 3);

if (size(sig,1) < size(sig,2)),
    warning(['Signal appears to be along rows, not columns.  Probably ' ...
             'should be transposed']);
end;

nonan = all(isfinite(sig),2);
if (~isempty(t))
    nonan = nonan & isfinite(t);
    t = t(nonan);
    if (any(abs((diff(t)-dt)/dt) > 1e-5)),
        error('Cannot handle unevenly sampling in time.');
    end;
end;
sig = sig(nonan,:);

N = size(sig,1);

y = fft(sig);
Py = y .* conj(y) / N;

f = 1/dt * (0:ceil(N/2))'/N;
Py = Py(1:length(f),:);

switch lower(opt.scale)
    case 'linear'
        % do nothing
    case 'db'
        Py = 10*log10(Py);
end

plot(f, Py);
xlabel('Frequency');
ylabel('Power');

if (nargout >= 2),
    varargout = {f,Py,y};
end;

