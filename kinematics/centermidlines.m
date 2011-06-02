function [x2,y2,hxctr,hyctr,ang] = centermidlines(x1,y1, indpeak)

nbt = size(indpeak,2);
nfr = size(x1,2);
npts = size(x1,1);

xctr0 = nans(2,nbt-1);
yctr0 = nans(2,nbt-1);
fr = zeros(1,nbt-1);

%run through each half tail beat
for ibt = 1:nbt-1,
    %get the center frame
    fr(ibt) = round((indpeak(end,ibt)+indpeak(end,ibt+1))/2);
    
    %and the center position of the head and tail
    xctr0(:,ibt) = (x1([1 end],indpeak(end,ibt)) + x1([1 end],indpeak(end,ibt+1)))/2;
    yctr0(:,ibt) = (y1([1 end],indpeak(end,ibt)) + y1([1 end],indpeak(end,ibt+1)))/2;
end;

%extrapolate linearly to the beginning and end of the defined midlines
span = first(all(isfinite(x1))):last(all(isfinite(x1)));

dxctr0 = diff(xctr0,[],2) ./ repmat(diff(fr),[2 1]);
dyctr0 = diff(yctr0,[],2) ./ repmat(diff(fr),[2 1]);

xstart = xctr0(:,1) - dxctr0(:,1).*((fr(1)-span(1))*[1;1]);
xend = xctr0(:,end) + dxctr0(:,end).*((span(end)-fr(end))*[1;1]);
xctr0 = [xstart xctr0 xend];
ystart = yctr0(:,1) - dyctr0(:,1).*((fr(1)-span(1))*[1;1]);
yend = yctr0(:,end) + dyctr0(:,end).*((span(end)-fr(end))*[1;1]);
yctr0 = [ystart yctr0 yend];
fr = [span(1) fr span(end)];

%now spline a smooth center position
xctr = nans(2,nfr);
yctr = nans(2,nfr);
xctr(1,span) = spline(fr,xctr0(1,:), span);
xctr(2,span) = spline(fr,xctr0(2,:), span);
yctr(1,span) = spline(fr,yctr0(1,:), span);
yctr(2,span) = spline(fr,yctr0(2,:), span);

%center each on the mean head position
hxctr = xctr(1,:);
hyctr = yctr(1,:);
x1ctr = x1 - repmat(hxctr,[npts 1]);
y1ctr = y1 - repmat(hyctr,[npts 1]);

%estimate midline angle
ang = atan2(diff(yctr), diff(xctr));

%rotate each midline to horizontal
ang1 = repmat(ang,[size(x1,1) 1]);
x2 =  x1ctr.*cos(ang1) + y1ctr.*sin(ang1);
y2 = -x1ctr.*sin(ang1) + y1ctr.*cos(ang1);

