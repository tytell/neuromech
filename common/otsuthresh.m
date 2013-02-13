function [B,t] = otsuthresh(I,block)
% function t = otsuthresh(I,blocksize)
% Calculates optimal threshold for an image (potentially separated into
% distinct blocks with size blocksize) according to the procedure from
% Otsu, N. (1978). Discriminant and least squares threshold selection. 
% In Proceedings of the Fourth International Joint Conference on Pattern
% Recognition, pp. 592-596.
%
% B is the thresholded image and t are the thresholds.
%
% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2010, Eric Tytell <tytell at jhu dot edu>


I = im2double(I);

if (nargin == 2),
	C = im2col(I,block,'distinct');
	N = prod(block);
else
	C = I(:);
	N = numel(I);
end;

L = 256;

[H,x] = hist(C,L);
H = shiftdim(H);
H = H / N;					% normalized histogram
i = linspace(0,1,L)';

omega = cumsum(H);
mu = cumsum(repmat(i,[1 size(H,2)]).*H);

muT = mu(L,:);
sigmaB = (repmat(muT,[L 1]) .* omega - mu).^2./(omega .* (1-omega));

[~,t] = max(sigmaB);
t = x(round(t));

C = C > repmat(t,[N 1]);

if (nargin == 2),
	B = col2im(C,block,size(I),'distinct');
else
	B = reshape(C,size(I));
end;



