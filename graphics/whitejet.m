function J = whitejet(m)
% J = whitejet(m)
% Like the jet colormap, but with white in the middle instead of green.

if nargin < 1
   m = size(get(gcf,'colormap'),1);
end

ind = ([1; 8; 24; 29; 35; 40; 56; 64]-1)/63;
col = [0 0 143; 0 16 255; 0 255 255; 255 255 255; ...
    255 255 255; 255 255 0; 239 0 0; 128 0 0] / 255;

J = zeros(m,3);
J(:,1) = interp1(round(ind*m)+1,col(:,1), (1:m)');
J(:,2) = interp1(round(ind*m)+1,col(:,2), (1:m)');
J(:,3) = interp1(round(ind*m)+1,col(:,3), (1:m)');

