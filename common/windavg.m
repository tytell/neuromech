function varargout = windavg(t,y,wind,twind,fcn)
% function ywind = windavg(t,y,wind,twind,fcn)
%    or    [ywind,twind] = windavg(...)

if ((nargin < 4) || isempty(twind)),
    tstart = floor(min(t(:))/wind)*wind - binsize/2;
    tend = ceil(max(t(:))/wind)*wind + binsize/2;

    twind = tstart:wind:tend;
end;

nt = length(t);

ind = zeros(0,length(twind));

a = 1;
b = 1;
for i = 1:length(twind),
    while ((a <= nt) && t(a) < twind(i)-wind),
        a = a+1;
    end;
    while ((b <= nt) && (t(b) < twind(i)+wind)),
        b = b+1;
    end;
    
    ind(1:b-a,i) = (a:b-1)';
end;

good = ind ~= 0;

yy = NaN(size(ind));
yy(good) = y(ind(good));

if (nargin < 5),
    fcn = @nanmean;
end;

ywind = feval(fcn,yy);

if (nargout == 1),
    varargout = {ywind};
else
    varargout = {ywind,twind};
end;


