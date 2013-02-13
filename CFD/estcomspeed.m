function [width,area,sn, comx,comy,comvelx,comvely] = ...
    estcomspeed(t,xms,yms,xns,yns,xls,xrs,yls,yrs)
% function [width,area,sn, comx,comy,comvelx,comvely] = ...
%    estcomspeed(t,xms,yms,xns,yns,xls,xrs,yls,yrs)
% Estimates the position and velocity of the center of mass, given the four
% sets of points that define the swimmer.

npt = size(xms,1);
nfr = size(xms,2);

width = sqrt((xls - xrs).^2 + (yls - yrs).^2);

sn = repmat(sqrt((xns(1,:)-xms(1,:)).^2 + (yns(1,:)-yms(1,:)).^2),[npt 1]) + ...
    [zeros(1,nfr); cumsum(sqrt(diff(xns).^2 + diff(yns).^2))];
len = sum(sqrt(diff(xms).^2 + diff(yms).^2));

sn = [zeros(1,nfr); sn(1:end-1,:); len];

%integrate to get total area
area = zeros(1,nfr);
for i = 1:nfr,
    area(i) = trapz(sn(:,i),[0; width(1:end-1,i); 0]);
end;

%integrate(x(s) * width(s) ds)/area
% to get COM position
comx = zeros(1,nfr);
comy = zeros(1,nfr);
for i = 1:nfr,
    comx(i) = trapz(sn(:,i),...
        [0; xns(1:end-1,i).*width(1:end-1,i); 0] ./ area(i));
    comy(i) = trapz(sn(:,i),...
        [0; yns(1:end-1,i).*width(1:end-1,i); 0] ./ area(i));
end;

%then derivative to get velocity
comvelx = deriv(t,comx);
comvely = deriv(t,comy);
