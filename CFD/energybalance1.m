function [totalke1,dissip1,boundflux1] = energybalance1(V0,opt)
% function [totalke1,dissip1,boundflux1] = energybalance1(V0,opt)
% Estimates total kinetic energy, dissipation and boundary flux in a fluid
% domain.  V0 comes from importsamrai.
%
% Units - totalke1 is [g/cm^3] * ([cm/s])^2 * (cm)^2 = g cm / s^2 = 
% ergs / cm
% dissip1 is [g / cm / s] * ([cm/s]/[cm])^2 * (cm)^2 = g cm / s^3 = 
% ergs / cm / s
% boundflux1 is [cm] * [g/cm^3] * ([cm/s])^2 * [cm/s] = g cm / s^3 =
% ergs / cm / s

totalke1 = 0.5 * opt.rho * trapzsamrai(V0, @(x) (x.U_0.^2 + x.U_1.^2));
dissip1 = 0.5 * opt.mu * trapzsamrai(V0, ...
    @(x) (4*x.dudx.^2 + 2*(x.dudy + x.dvdx).^2 + 4*x.dvdy.^2));

lo = min(cat(2,V0.xlo),[],2);
hi = max(cat(2,V0.xup),[],2);
dxhi = max(cat(2,V0.dx),[],2);
dxlo = min(cat(2,V0.dx),[],2);
dyhi = max(cat(2,V0.dy),[],2);
dylo = min(cat(2,V0.dy),[],2);

boundy0 = lo(2)+dyhi:dylo:hi(2)-dyhi;
boundy0 = boundy0 + (hi(2)-dyhi - boundy0(end))/2;
boundx0 = lo(1)+dxhi:dxlo:hi(1)-dxhi;
boundx0 = boundx0 + (hi(1)-dxhi - boundx0(end))/2;

boundx1 = NaN(4,max(length(boundx0),length(boundy0)));
boundy1 = NaN(size(boundx1));
boundx1(1,1:length(boundy0)) = lo(1)+dxhi;
boundy1(1,1:length(boundy0)) = boundy0;
boundx1(2,1:length(boundx0)) = boundx0;
boundy1(2,1:length(boundx0)) = lo(2)+dyhi;
boundx1(3,1:length(boundy0)) = hi(1)-dxhi;
boundy1(3,1:length(boundy0)) = boundy0;
boundx1(4,1:length(boundx0)) = boundx0;
boundy1(4,1:length(boundx0)) = hi(2)-dyhi;

boundvals = interpsamrai(V0, boundx1,boundy1, 'vars',{'U_0','U_1'});
boundke = 0.5 * opt.rho * (boundvals.U_0.^2 + boundvals.U_1.^2);
good = ~isnan(boundx1);
boundflux1(1) = trapz(boundy1(1,good(1,:)),...
    boundke(1,good(1,:)).*boundvals.U_0(1,good(1,:)));
boundflux1(2) = -trapz(boundx1(2,good(2,:)),...
    boundke(2,good(2,:)).*boundvals.U_1(2,good(2,:)));
boundflux1(3) = -trapz(boundy1(3,good(3,:)),...
    boundke(3,good(3,:)).*boundvals.U_0(3,good(3,:)));
boundflux1(4) = trapz(boundx1(4,good(4,:)),...
    boundke(4,good(4,:)).*boundvals.U_1(4,good(4,:)));

