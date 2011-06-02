function varargout = plotfftpower(varargin)

if ((nargin >= 2) & isnumeric(varargin{1}) & isnumeric(varargin{2})),
    if (length(varargin{2}) == 1),
        dt = 1/varargin{2};
        sig = varargin{1};
    else,
        t = varargin{1};
        sig = varargin{2};

        dt = t(2)-t(1);
    end;
elseif (nargin == 1),
    sig = varargin{1};
    dt = 1;
end;

if (size(sig,1) < size(sig,2)),
    warning(['Signal appears to be along rows, not columns.  Probably ' ...
             'should be transposed']);
end;

nonan = all(isfinite(sig),2) & isfinite(t);
t = t(nonan);
if (any(abs((diff(t)-dt)/dt) > 1e-5)),
    error('Cannot handle unevenly sampling in time.');
end;
sig = sig(nonan,:);

N = size(sig,1);

y = fft(sig);
Py = y .* conj(y) / N;

f = 1/dt * (0:ceil(N/2))'/N;
Py = Py(1:length(f),:);

plot(f, Py);
xlabel('Frequency');
ylabel('Power');

if (nargout == 2),
    varargout = {f,Py};
end;

