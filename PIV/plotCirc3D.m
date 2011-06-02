function [h,hbody,maxcirc] = rgnIsoPlot(vxr,pivframes,bts, t,flowvel,...
                                sheetpos,lenmm, ...
                                ctrx,ctry, bodyx,bodyy)
% Plots vortices as surfaces through time, coloring by circulation
% Works using the vxr output of identifyVortices

if (nargin == 7),
    ctrx = zeros(1,size(sheetpos,2));
    ctry = zeros(1,size(sheetpos,2));
    bodyx = [];
    bodyy = [];
elseif (isempty(ctrx)),
    ctrx = zeros(1,size(bodyx,2));
    ctry = zeros(1,size(bodyx,2));
end;

if (all(ctrx == 0)),
    iscentered = false;
else
    iscentered = true;
end;

if (~ishold),
    cla;
end;
hold on;

% number of beats to work from
nphase = size(pivframes,1);
npiv = sum(flatten(isfinite(pivframes)));

pfr = reshape(pivframes,[nphase/2 size(pivframes,2)*2]);
bts2 = flatten(([bts; bts+0.5]-1)*2)' + 1;

btframes = flatten(pfr(:,bts2))';
btpiv = sub2ind(size(pfr),repmat((1:nphase/2)',[1 length(bts2)]), ...
                repmat(bts2,[nphase/2 1]));
btpiv = flatten(btpiv)' - sum(isnan(pivframes(:,1)));

k = find((btpiv >= 1) & (btpiv <= npiv));
btpiv = btpiv(k);
btframes = btframes(k);

phi = linspace(0,2*pi,17)';

% make the artificial z position, including any known changes in position
% relative to the sheet.  Start upstream by the first sheet pos - when
% zpos reaches 0, then stuff at the first sheet pos should have convected
% to the tail
zpos = repmat(NaN,[1 npiv]);
zpos(btpiv) = cumtrapz(t(btframes), flowvel + ...
                       deriv(t(btframes),sheetpos(btframes))) ...
    + lenmm - sheetpos(btframes(1));

h = -1*ones(size(vxr));
maxcirc = NaN*ones(size(vxr));
for ii = 1:length(vxr),
    i = ii;

    k = find((vxr(i).frames >= btpiv(1)) & ...
             (vxr(i).frames <= btpiv(end)));

    vxrbtfr = vxr(i).frames(k) - btpiv(1) + 1;
    vxrpivfr = vxr(i).frames(k);

    %adjust angle to vary as smoothly as possible, mod pi (not 2*pi)
    %first find the mean angle
    cc = cos(vxr(i).angle(k));
    ss = sin(vxr(i).angle(k));
    meanang = atan2(nanmean(ss),nanmean(cc));
    
    ang1 = mod(vxr(i).angle(k) - meanang,pi) + meanang;

    a = vxr(i).majordiam(k)/2;
    b = vxr(i).minordiam(k)/2;
    ang2 = atan2(a.*sin(ang1), ...
                 b.*cos(ang1));

    x1 = zeros(length(phi),length(k));
    y1 = zeros(length(phi),length(k));
    for j = 1:length(k),
        x1(:,j) = a(j)*cos(phi-ang2(j)).*cos(ang1(j)) - ...
                  b(j)*sin(phi-ang2(j)).*sin(ang1(j)) + ...
                  vxr(i).ctrx(j);
        y1(:,j) = a(j)*cos(phi-ang2(j)).*sin(ang1(j)) + ...
                  b(j)*sin(phi-ang2(j)).*cos(ang1(j)) + ...
                  vxr(i).ctry(j);
    end;
    
    
    z1 = repmat(zpos(vxrpivfr),[size(x1,1) 1]);

    if ((length(k) > 1) & ~isempty(x1)),
        x1 = x1 - repmat(ctrx(btframes(vxrbtfr)), ...
                         [size(x1,1) 1]);
        y1 = y1 - repmat(ctry(btframes(vxrbtfr)), ...
                         [size(y1,1) 1]);

        h(i) = surf(x1,y1,z1,repmat(vxr(i).circmax(k),[size(x1,1) 1]), ...
                    'EdgeColor','none');
        maxcirc(i) = max(vxr(i).circmax(k));
    end;
end;

if (~isempty(bodyx)),
    hbody = surf(bodyx(:,btframes)' - ...
                 repmat(ctrx(btframes),[size(bodyx,1) 1])', ...
                 bodyy(:,btframes)' - ...
                 repmat(ctry(btframes),[size(bodyy,1) 1])', ...
                 repmat(zpos(btpiv),[size(bodyx,1) 1])', ...
                 'FaceColor','k','EdgeColor','k');

    alpha(hbody,0.5);
else
    hbody = -1;
end;

axis equal vis3d tight;

clear x1 y1 z1 cc ss meanang ang1 ang2 phi a b i k j bts btframes zpos ...
      good ii;

hold off;



