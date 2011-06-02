function varargout = spiketrigger(spikeind, rng, A, varargin)

p = 1;
while ((p <= length(varargin)) && isnumeric(varargin{p}) && ...
        (ndims(varargin{p}) == ndims(A)) && all(size(varargin{p}) == size(A)))
    p = p+1;
end;

C = [A varargin(1:p-1)];
N = length(C);

off = (rng(1):rng(2))';
noff = length(off);

spikeind = makerow(spikeind);

ind = spikeind(ones(noff,1),:) + off(:,ones(1,length(spikeind)));
good = all((ind >= 1) & (ind <= length(A)));

for i = 1:N,
    A1 = NaN(size(ind));
    A1(:,good) = C{i}(ind(:,good));
    
    C{i} = A1;
end;

varargout = C;
