function [B,t] = otsuthresh(I,block)
% function t = otsuthresh(I)

I = im2double(I);

if (nargin == 2),
	C = im2col(I,block,'distinct');
	N = prod(block);
else
	C = I(:);
	N = prod(size(I));
end;

L = 256;

H = hist(C,L);
H = shiftdim(H);
H = H / N;					% normalized histogram
i = linspace(0,1,L)';

omega = cumsum(H);
mu = cumsum(repmat(i,[1 size(H,2)]).*H);

muT = mu(L,:);
sigmaB = (repmat(muT,[L 1]) .* omega - mu).^2./(omega .* (1-omega));

[s,t] = max(sigmaB);
t = t/L;

C = C > repmat(t,[N 1]);

if (nargin == 2),
	B = col2im(C,block,size(I),'distinct');
else
	B = reshape(C,size(I));
end;



