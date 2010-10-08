function f = unmod(f,m,delta,dim)
% function y = unmod(f,m,delta)
% Removes the modulus from discrete data, by looking for large jumps and
% assuming they represent points when the values have wrapped around in
% modulus.  Similar to unwrap, but can work with any modulus.
%
% delta is optional and represents the size of the jump.  Default is m/10
%
% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2010, Eric Tytell <tytell at jhu dot edu>

% deal with optional arguments
if (nargin < 4)
	dim = [];
	if (nargin < 3),
		delta = [];
		if (nargin < 2),
			m = 2*pi;
		end;       
	end;
end;

if (isempty(dim)),
    if ((ndims(f) == 2) && (size(f,1) == 1)),
        dim = 2;
    else
        dim = 1;
    end;
end;

if (isempty(f) || (size(f,dim) < 2)),
    return;
end;

if (isempty(m)),
    m = 2*pi;
end;

if (isempty(delta))
	delta = m/10;
end;

sz = size(f);
otherdims = [1:dim-1 dim+1:length(sz)];
rest = prod(sz(otherdims));

f = permute(f, [dim otherdims]);
f = reshape(f, [sz(dim) rest]);

% step through any dimensions of f beyond the first
for i = 1:rest,
	% only operate on the finite elements
	good = isfinite(f(:,i));
	
	% take their differences
    if (sum(good) > 1),
        f1 = f(good,i);
        d = diff(f1);

        % eliminate differences that are close to m
        jump = (abs(m - abs(d)) < delta) & (d > 0);
        d(jump) = d(jump) - m;
        jump = (abs(m - abs(d)) < delta) & (d < 0);
        d(jump) = d(jump) + m;

        % use a cumulative sum to get f back
        f(good,i) = cumsum([f1(1); d]);
    end;
end;

f = reshape(f,[sz(dim) sz(otherdims)]);
f = permute(f,[2:dim 1 dim+1:length(sz)]);
