function [A,filt] = smooth2d(A, varargin)
% function [A,filt] = smooth2d(A, ord,c1,c2)
%   or                smooth2d(A, filt)
% Simple 2d smoothing based on a Hamming window.  ord is the filter order,
% c1 and c2 are the cutoffs along dimension 1 and 2, respectively.
%
% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2010, Eric Tytell

if (nargin == 4)
    [ord, c1,c2] = varargin{:};

    [f1,f2] = freqspace(ord,'meshgrid');
    
    Hd = zeros(ord,ord);
    d = sqrt((f1/c1).^2 + (f2/c2).^2) <= 1;
    Hd(d) = 1;
    
    filt = fwind1(Hd,hamming(ord));
elseif (nargin == 2)
    [filt] = varargin{:};
else
    error('smooth2d:args','Wrong number of arguments');
end;

A = filter2(filt,A);
