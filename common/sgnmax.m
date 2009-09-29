function [m,ind] = sgnmax(varargin)

[hi,indhi] = max(varargin{:});
[lo,indlo] = min(varargin{:});

ishi = abs(hi) > abs(lo);
m = lo;
ind = indlo;
m(ishi) = hi;
ind(ishi) = indhi;

