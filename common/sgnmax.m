function [m,ind] = sgnmax(varargin)
% function [m,ind] = sgnmax(varargin)
% Returns the values with the maximum magnitude/absolute value, but with
% the signs correct.
%
% sgnmax([5 -5 3; 2 3 -8]) = [5 -5 -8]
%
% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2010, Eric Tytell <tytell at jhu dot edu>


[hi,indhi] = max(varargin{:});
[lo,indlo] = min(varargin{:});

ishi = abs(hi) > abs(lo);
m = lo;
ind = indlo;
m(ishi) = hi;
ind(ishi) = indhi;

