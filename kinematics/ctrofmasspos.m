function [comx,comy] = ctrofmasspos(mx,my, width,height)

if ((nargin == 2) || isempty(width)),
    %just the outline
    %calculate midline and width
    ox = mx;
    oy = my;
    
    npts = size(ox,1)/2;
    
    left = 1:npts;
    right = npts*2:-1:npts+1;
    
    mx = (ox(left,:) + ox(right,:))/2;
    my = (oy(left,:) + oy(right,:))/2;
    
    width = sqrt((ox(left,:)-ox(right,:)).^2 + (oy(left,:)-oy(right,:)).^2)/2;
    width = nanmedian(width,2);
end;

sz = size(mx);
mx = reshape(mx,[sz(1) prod(sz(2:end))]);
my = reshape(my,[sz(1) prod(sz(2:end))]);

npts = size(mx,1);
nfr = size(mx,2);

if ((nargin < 4) || isempty(height)),
    %no height -- assume constant height
    height = ones(npts,1) / npts;
end;

width = makecol(width);
height = makecol(height);

%normalize width and height to sum to 1
whsum = sum(width .* height);
width = repmat(width, [1 nfr]);
height = repmat(height, [1 nfr]);

%center of mass is just mean position, weighted by width and height
comx = sum(mx .* width .* height / whsum);
comy = sum(my .* width .* height / whsum);

comx = reshape(comx,[1 sz(2:end)]);
comy = reshape(comy,[1 sz(2:end)]);


   