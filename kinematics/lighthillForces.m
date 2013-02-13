function [Fthrust,Fthrustmean0,Flateral, Pwaste] = lighthillForces(massperlen,s,t,swimvel,x,y)
% function [Fthrust,Fthrustmean0,Flateral] = 
%                lighthillForces(massperlen,s,t,swimvel,dxdt,dydt,dxds,dyds)
%
% Calculates swimming forces according to the Lighthill large amplitude
% elongated body theory (1971).  Note that signs are very important.  Here,
% we have x and s the inverse of Lighthill (he uses a, the distance from the
% tail, where we use s, the distance from the head).  Our fish is swimming
% in the negative x direction.  But, pass swimvel in as a positive number.
% Our y is identical to his z.
%
% The equations are (in Lighthill's coordinates)
%	Fthrust = (m w dz/dt - 0.5 m w^2 dx/da)_(a=0) - 
%			d/dt ( Integrate(-m w dz/da  da, a = 0 to l) )
%	Flateral = (-m w dx/dt - 0.5 m w^2 dz/da)_(a=0) -
%			d/dt ( Integrate(m w dx/da  da, a = 0 to l) )
%
% where w is the velocity perpedicular to the midline,
%	w = dz/dt dx/da - dx/dt dz/da
%
% We convert the data into Lighthill's coordinate system, evaluate the
% equations there, and then convert back.  Note that thrust should be in the
% same direction as the swimming velocity.
%
% Units:
%   velperp -> [cm/s] * [cm / cm] -> cm/s
%   integ = [cm] * [g/cm] * [cm/s] * [cm/cm] = g cm / s
%   Flateral = [g/cm] * [cm/s] * (-[cm/s] - [cm/s]*[cm/cm]) + d integ / dt
%            = [g cm / s^2] + d [g cm/s] / dt = [g cm / s^2] = dynes

% Convert from my coordinate system to Lighthill's
if (range(x(:)) < range(y(:))),
    ishoriz = 0;
    swap = -y;
    y = x;
    x = swap;
else
    ishoriz = 1;
end;

a = s(end) - s;
x = -x;

dyda = deriv(a,y);
% NB: dxda > 0, because we're swimming in the positive x direction (now
% that we flipped the sign on x) and a increases from tail to head, so if
% the head is pointing in the direction we're swimming, then dxda > 0.
% By definition (dx/da)^2 + (dy/da)^2 = 1.  This means the midline always
% has the same length.
dxda = sqrt(1 - dyda.^2);

if (size(swimvel,2) == 1)
    swimvel = repmat(swimvel,[1 size(x,2)]);
end;

dxdt = deriv(t,x, 2) + swimvel(ones(size(x,1),1),:);
dydt = deriv(t,y, 2);

% velocity tangent to the midline
velperp = dydt.*dxda - dxdt.*dyda;

Fthrust = massperlen(end).*velperp(end,:) .* ...
          (dydt(end,:) - 0.5*velperp(end,:).*dxda(end,:));

% note a(1) = head and a(end) = tail.  We are supposed to integrate
% from tail to head, which leads to the negative.
integ = -trapz(a,massperlen(:,ones(1,size(velperp,2))).*velperp.*dyda);

Fthrustmean0 = deriv(t,integ);

% first term of the lateral force
Flateral = massperlen(end).*velperp(end,:) .* ...
    (-dxdt(end,:) - 0.5*velperp(end,:).*dyda(end,:));

% second term
integ = trapz(a,massperlen(:,ones(1,size(velperp,2))).*velperp.*dxda);
Flateral = Flateral + deriv(t,integ);

% calculate wasted power
veltan = dxdt.*dxda + dydt.*dyda;
Pwaste = 0.5*massperlen(end).*velperp(end,:).^2.*veltan(end,:);

Pbody = 0.5*massperlen(:,ones(1,size(velperp,2))).*velperp.^2.*veltan;

% convert coordinate systems again
Fthrust = -Fthrust;
Fthrustmean0 = -Fthrustmean0;

% dxdt = dxdt - swimvel;
% veltan = dydt.*dxds - dxdt.*dyds;
% 
% Fthrust = massperlen.*veltan(end,:) .* ...
% 			(dydt(end,:) - 0.5*veltan(end,:).*dxds(end,:));
% 
% Ft = trapz(s,massperlen.*veltan.*dyds);
% Fthrustmean0 = -deriv(t,Ft);
% 
% Flateral = massperlen.*veltan(end,:) .* ...
% 			(dxdt(end,:) + 0.5*veltan(end,:).*dyds(end,:));
% 
% Ft = trapz(s,massperlen.*veltan.*dxds);
% Flateral = Flateral - deriv(t,Ft);
