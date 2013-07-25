function varargout = angmeantest(varargin)
% ANGMEANTEST - Rayleigh test for a significant mean angle
%
% function [P,z] = angmeantest(ang)
%   or             angmeantest(angx,angy)
% 
% where ang is a vector of angles in radians, or (angx,angy) are the
% coordinates of a set of unit vectors.
% With no output arguments, it displays a table output.

% Mercurial revision hash: $Revision$ $Date$
% See http://rcn.ccs.tulane.edu/index.php5/Tytell_Matlab
% Copyright (c) 2012, Eric Tytell <tytell at jhu dot edu>

if (nargin == 1),
  ang = varargin{1};

  angx = nanmean(cos(ang));
  angy = nanmean(sin(ang));
  n = sum(isfinite(ang));
elseif (nargin == 2),
  angx = nanmean(varargin{1});
  angy = nanmean(varargin{2});
  n = sum(isfinite(varargin{1}));
end;

a = atan2(angy,angx);
r = sqrt(angx.^2 + angy.^2);
z = n.*r.^2;
s = sqrt(-2*log(r));

P = raylzinv(z,n);

if (nargout == 0),
    fprintf('Rayleigh test against random orientation.\n');
    fprintf('%8s %8s %8s %8s %8s %8s %8s\n','Sample','C','S','ang','R','z','P');
    fprintf('%8d %8.3f %8.3f %8.3f %8.3f %8.3f %8.5f\n',[1:size(r,2); angx; angy; ...
        a*180/pi; r; z; P]);
else
    varargout = {P,z,a,r,s};
end;
