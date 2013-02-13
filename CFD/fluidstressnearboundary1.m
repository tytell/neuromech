function S1 = fluidstressnearboundary1(V, boundx0,boundy0,boundu0,boundv0, opt)
% function S1 = fluidstressnearboundary1(V, boundx0,boundy0,boundu0,boundv0, opt)
% Given the boundary boundx0,boundy0,boundu0,boundv0, and the adaptive grid
% data in V (representing a single time point), estimate the jump in fluid
% stress across the boundary, using the "fluid stress" method of Williams
% et al. 2009.

npt = 321;

levelnum = cat(2,V.level_number);
isfinest = levelnum == max(levelnum);
V = V(isfinest);

d = V(1).dx;

s0 = [0; cumsum(sqrt(diff(boundx0).^2 + diff(boundy0).^2))];
stail = s0(npt);

%smallest two values that are greater than d and divide each side evenly
d1 = stail / floor(stail/d);
d2 = (s0(end)-stail) / floor((s0(end)-stail)/d);

s1 = [0:d1:stail stail+d2:d2:s0(end)];

boundx = interp1(s0,boundx0, s1, 'spline');
boundy = interp1(s0,boundy0, s1, 'spline');
boundux = interp1(s0,boundu0, s1, 'spline');
bounduy = interp1(s0,boundv0, s1, 'spline');

tanx1 = zeros(size(boundx));
tany1 = zeros(size(boundy));

tanx1(2:end-1) = boundx(3:end) - boundx(1:end-2);
tany1(2:end-1) = boundy(3:end) - boundy(1:end-2);
tanx1([1 end]) = boundx(2) - boundx(end-1);
tany1([1 end]) = boundy(2) - boundy(end-1);

mag = sqrt(tanx1.^2 + tany1.^2);
tanx1 = tanx1 ./ mag;
tany1 = tany1 ./ mag;

normx1 = tany1;
normy1 = -tanx1;

S1.tanx = tanx1;
S1.tany = tany1;

boundun = boundux.*normx1 + bounduy.*normy1;
boundut = boundux.*tanx1 + bounduy.*tany1;

n1 = d*(0:5);

[ss,nn] = meshgrid(s1,n1);
S1.s = ss;
S1.n = nn;

nearboundx1 = repmat(boundx, [length(n1) 1]) + repmat(normx1, [length(n1) 1]) .* nn;
nearboundy1 = repmat(boundy, [length(n1) 1]) + repmat(normy1, [length(n1) 1]) .* nn;
S1.nearboundx = nearboundx1;
S1.nearboundy = nearboundy1;

%distance between neighboring points along the tangential direction is
%different depending on where you are on the normal axis
tdist = sqrt(diff(nearboundx1,[],2).^2 + diff(nearboundy1,[],2).^2);

tanx2 = repmat(tanx1, [length(n1) 1]);
tany2 = repmat(tany1, [length(n1),1]);
normx2 = repmat(normx1, [length(n1) 1]);
normy2 = repmat(normy1, [length(n1),1]);

xlo = cat(2,V.xlo);
xup = cat(2,V.xup);

ylo = xlo(2,:);
yup = xup(2,:);
xlo = xlo(1,:);
xup = xup(1,:);

nearboundut = NaN(size(nearboundx1));
nearboundun = NaN(size(nearboundy1));
nearboundp = NaN(size(nearboundx1));

% estimate the tangential and normal velocities close to the boundary and
% the pressure near the boundary
progressval = 0;
for i = 1:length(V),
    inrng = (nearboundx1 >= xlo(i)) & (nearboundx1 <= xup(i)) & ...
        (nearboundy1 >= ylo(i)) & (nearboundy1 <= yup(i));
    
    [xgrid,ygrid] = meshgrid(V(i).px,V(i).py);
    
    ux = V(i).U_0;
    uy = V(i).U_1;
    p = V(i).P;
    
    in = inpolygon(xgrid,ygrid, boundx0,boundy0);
    
    ux(in) = NaN;
    uy(in) = NaN;
    p(in) = NaN;
    
    pcubic = NaN(size(nearboundx1));
    plinear = NaN(size(nearboundx1));
    pcubic(inrng) = interp2(xgrid,ygrid, p, nearboundx1(inrng),nearboundy1(inrng), '*cubic');
    plinear(inrng) = interp2(xgrid,ygrid, p, nearboundx1(inrng),nearboundy1(inrng), '*linear');
    nearboundp(inrng & isfinite(plinear)) = plinear(inrng & isfinite(plinear));
    nearboundp(inrng & isfinite(pcubic)) = pcubic(inrng & isfinite(pcubic));
    
    ux1 = NaN(size(nearboundx1));
    uy1 = NaN(size(nearboundx1));
    
    ux1(inrng) = interp2(xgrid,ygrid, ux, nearboundx1(inrng),nearboundy1(inrng), '*linear');
    uy1(inrng) = interp2(xgrid,ygrid, uy, nearboundx1(inrng),nearboundy1(inrng), '*linear');
    
    ut1 = ux1.*tanx2 + uy1.*tany2;
    un1 = ux1.*normx2 + uy1.*normy2;
    
    nearboundut(inrng) = ut1(inrng);
    nearboundun(inrng) = un1(inrng);
    nearboundut(1,inrng(1,:)) = boundut(inrng(1,:));
    nearboundun(1,inrng(1,:)) = boundun(inrng(1,:));
    
    if (floor(i/length(V)/0.1) ~= progressval)
        p = floor(i/length(V)/0.1);
        fprintf('%d%% ',p*10);
        progressval = p;
    end;
end;
fprintf('\n');

for i = 1:length(n1),
    good = ~isnan(nearboundun(i,:));
    nearboundun(i,~good) = interp1(s1(good), nearboundun(i,good), s1(~good), 'cubic');
    nearboundut(i,~good) = interp1(s1(good), nearboundut(i,good), s1(~good), 'cubic');
    
    good = ~isnan(nearboundp(i,:));
    if (sum(good) > 0.5*length(good))
        nearboundp(i,~good) = interp1(s1(good), nearboundp(i,good), s1(~good), 'cubic');
    end;
end;

dundn = deriv(nn,nearboundun);
dutdn = deriv(nn,nearboundut);

%take the derivative along the tangential direction, correctly accounting
%for stretching as one moves away on the normal axis.  Also take into
%account that the contour wraps around
dundt = NaN(size(nearboundun));
dundt(:,2:end-1) = (nearboundun(:,3:end)-nearboundun(:,1:end-2)) ./ ...
    (tdist(:,1:end-1) + tdist(:,2:end));
dundt(:,1) = (nearboundun(:,2)-nearboundun(:,end-1)) ./ ...
    (tdist(:,1) + tdist(:,end));
dundt(:,end) = (nearboundun(:,2)-nearboundun(:,end-1)) ./ ...
    (tdist(:,1) + tdist(:,end));

S1.nearboundp = nearboundp;
S1.nearboundun = nearboundun;
S1.nearboundut = nearboundut;
S1.vorticity = dundt - dutdn;

% mu [g/cm/s] * dundn [cm/s/cm] -> g/cm/s^2 = dynes/cm^2
% p [dynes/cm^2] = [g cm /s^2 / cm^2] = [g / cm / s^2]
S1.normstress = 2*opt.mu*dundn - nearboundp;
S1.tanstress = opt.mu*(dutdn + dundt);
