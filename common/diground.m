function r = diground(a,dig)
% function r = diground(a,dig)
% Rounds to a specific digit.
%   diground(10.268,0.01) = 10.27
%   diground(10.268,0.1) = 10.3
%   diground(10.268,0.25) = 10.25
%
% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2010, Eric Tytell <tytell at jhu dot edu>

r = dig*round(a/dig);
