function theta = VMrand(mu,k,varargin)
% VMRAND - Random angles from a von Mises distribution
%
% function theta = VMrand(mu,k,m,n,...)
%   or     theta = VMrang(mu,k,[m n ...])
% 
% Generates random angles ([0 2pi]) according to a von Mises distribution with mean
% mu and concentration k.  Size of the result is [m n ...]

% Mercurial revision hash: $Revision$ $Date$
% See http://rcn.ccs.tulane.edu/index.php5/Tytell_Matlab
% Copyright (c) 2012, Eric Tytell <tytell at jhu dot edu>


switch length(varargin)
    case 0,
        sz = [1 1];
    case 1,
        sz = [varargin{1}(:)' 1];
    otherwise,
        sz = cat(1,varargin{:});
end;

N = prod(sz);

a = 1 + sqrt(1 + 4*k^2);
b = (a - sqrt(2*a))/(2*k);
r = (1 + b^2)/(2*b);

theta = zeros(sz);

for i = 1:N,
  done = 0;
  while ~done,
    U1 = rand;
    z = cos(pi*U1);
    f = (1 + r*z)/(r + z);
    c = k*(r-f);

    U2 = rand;
    if (c*(2-c) - U2 > 0),
      done = 1;
    elseif ((log(c/U2) + 1 - c) >= 0),
      done = 1;
    end;
  end;
  U3 = rand;
  theta(i) = mod(sign(U3 - 0.5)*acos(f) + mu, 2*pi);
end;

 