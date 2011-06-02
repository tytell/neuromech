function dcev = disccomplexeig(varargin)
% dcev = disccomplexeig(DU)
%      or disccomplexeig(x,y,u,v)
% Returns the discriminant for complex eigenvalues based on the velocity
% derivative matrix, DU.  The discriminant is
%		(du/dx + dv/dy)^2 - 4*(du/dx dv/dy - du/dy dv/dx)
% and thus its units are 1/time^2.
% To make things more similar to DaVis's "swirling strength", I flip the
% sign and get rid of anything below zero

if (nargin == 4),
    u = varargin{3};
    v = varargin{4};
    if (ndims(u) > 3),
        sz = size(u);
        u = flatten(u,3:ndims(u));
        v = flatten(v,3:ndims(u));
    end;
    DU = velderiv(varargin{1:2},u,v);
else
    DU = varargin{1};
end;

dcev = (cat(3,DU.dudx) + cat(3,DU.dvdy)).^2 - ...
       4*(cat(3,DU.dudx) .* cat(3,DU.dvdy) - ...
          cat(3,DU.dudy) .* cat(3,DU.dvdx));
		
dcev = -dcev;
dcev(dcev < 0) = 0;

if (exist('sz')),
    dcev = reshape(dcev,sz);
end;

