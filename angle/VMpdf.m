function p = VMpdf(theta,mu,k)
% function p = VMpdf(theta,mu,k)
% PDF for the von Mises distribution.  Returns probabilities p of angles
% theta, given a mean mu and concentration k.
%
% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2010, Eric Tytell

p = 1./(2*pi*besseli(0,k)) .* exp(k.*cos(theta - mu));
