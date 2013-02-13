function [good,replace,ur,vr] = validatePIV(u,v, medmax)
% function good = validatePIV(u,v, medmax)

filt = [-1   4  -6   4  -1; ...
         4 -16  24 -16   4; ...
        -6  24 100  24  -6; ...
         4 -16  24 -16   4; ...
        -1   4  -6   4  -1];
    
nbh = [5 5];
off = floor(nbh/2);
ctr = nbh(1)*(ceil(nbh(2)/2)-1) + ceil(nbh(1)/2);
sz = size(u);

ngoodneighbors = round(0.9 * numel(filt));

%turn u and v into columns containing the neighborhood around each point
ucol = zeros(nbh(1)*nbh(2), (sz(1)-nbh(1)+1)*(sz(2)-nbh(2)+1), sz(3));
vcol = zeros(nbh(1)*nbh(2), (sz(1)-nbh(1)+1)*(sz(2)-nbh(2)+1), sz(3));
for i = 1:size(u,3),
    ucol(:,:,i) = im2col(u(:,:,i),nbh, 'sliding');
    vcol(:,:,i) = im2col(v(:,:,i),nbh, 'sliding');
end;

%take the median without the center point
noctr = [1:ctr-1 ctr+1:prod(nbh)]';
umed = nanmedian(ucol(noctr,:,:));
vmed = nanmedian(vcol(noctr,:,:));

%no do a "standard deviation" = RMS around the median
urms = sqrt(nanmean((ucol(noctr,:,:) - repmat(umed,[length(noctr) 1 1])).^2));
vrms = sqrt(nanmean((vcol(noctr,:,:) - repmat(vmed,[length(noctr) 1 1])).^2));
urms(urms == 0) = NaN;
vrms(vrms == 0) = NaN;

%check how far the center is from the rest
udev = abs(ucol(ctr,:,:) - umed) ./ urms;
udev = reshape(udev, [sz(1)-nbh(1)+1 sz(2)-nbh(2)+1 sz(3)]);
vdev = abs(vcol(ctr,:,:) - vmed) ./ vrms;
vdev = reshape(vdev, [sz(1)-nbh(1)+1 sz(2)-nbh(2)+1 sz(3)]);

%keep ones that deviate less than medmax
good = true(size(u));
good(off(1)+1:end-off(1), off(2)+1:end-off(2), :) = (udev <= medmax) & (vdev <= medmax);

%now interpolate
goodcol = zeros(nbh(1)*nbh(2), (sz(1)-nbh(1)+1)*(sz(2)-nbh(2)+1), sz(3));
for i = 1:size(good,3),
    goodcol(:,:,i) = im2col(good(:,:,i),nbh, 'sliding');
end;

%replace if most of the neighbors are good
replace1 = ~goodcol(ctr,:,:) & (sum(goodcol(noctr,:,:)) >= ngoodneighbors);

%flatten out the frames for the moment
goodcol = flatten(goodcol,[2 3]);
ucol = flatten(ucol,[2 3]);
vcol = flatten(vcol,[2 3]);

%build up the filter matrix, but only for those that need replacing to save
%memory
filt = repmat(filt(:),[1 sum(replace1(:))]);

%don't apply the filter to bad points
filt(~goodcol(:,replace1) | ~isfinite(ucol(:,replace1))) = 0;
ucol(~isfinite(ucol)) = 0;
vcol(~isfinite(vcol)) = 0;

%should sum to 1
filtsum = repmat(sum(filt),[size(filt,1) 1]);
filt = filt ./ filtsum;

%apply the filter
ur = sum(ucol(:,replace1) .* filt);
vr = sum(vcol(:,replace1) .* filt);

%check to make sure the replacement values satisfy the median filter
goodreplace = (abs(ur - umed(:,replace1)) ./ urms(:,replace1) <= medmax) & ...
    (abs(vr - vmed(:,replace1)) ./ vrms(:,replace1) <= medmax);

replace1(replace1) = goodreplace;
ur = ur(goodreplace);
vr = vr(goodreplace);

%and reorganize the replacement matrix so that it's the same size as u
replace = false(size(u));
replace(off(1)+1:end-off(1), off(2)+1:end-off(2), :) = ...
    reshape(replace1, [sz(1)-nbh(1)+1 sz(2)-nbh(2)+1 sz(3)]);





    
    