function [x, jitter] = beeswarm(x, varargin)

opt.nbins = 20;
opt.edges = [];

opt = parsevarargin(opt, varargin, 2);

if ~isempty(opt.edges)
    [n,edges,bin] = histcounts(x, opt.edges);
else
    [n,edges,bin] = histcounts(x, opt.nbins);
end

jitter = zeros(size(x));
for i = 1:length(n)
    isbin = bin == i;
    jitter(isbin) = 0:sum(isbin)-1;
end

if nargout == 0
    plot(x, jitter, 'o');
end
