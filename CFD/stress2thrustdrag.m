function S = stress2thrustdrag(t,cycle,S, varargin)
% function S = stress2thrustdrag(t,cycle,S, varargin)
% Adds values to the structure S to convert the normal and tangential
% stresses into axial and lateral forces, then into "thrust" and "drag"
% impulse.  For simplicity, "thrust" is simply defined as positive in the
% axial direction, and "drag" is negative in the axial direction.
%
% Units:
%   S.Iaxialtot = integral S.Faxialtot dt = [dynes/cm] * [s] = [dynes*s/cm]
%   = impulse/cm

opt.optimizerotation = false;
opt.anglerange = [-20 20];
opt.steadycycle = [];

opt = parsevarargin(opt,varargin,3);

for fr = 1:length(t)
    s1 = S.s{fr}(1,:);
    S.Flateraltot(fr) = trapz(s1,S.lateralstresst{fr}+S.lateralstressn{fr});
end;

if (opt.optimizerotation && isempty(opt.steadycycle)),
    error('Must include steadycycle as an option for optimizing rotation');
end;

if (opt.optimizerotation)
    ang = (opt.anglerange(1):opt.anglerange(2)) * pi/180;
    dang = ang(2)-ang(1);
    
    Iax1 = zeros(cycle(end)-1, length(ang));
    ang0 = NaN(cycle(end)-1,1);
    
    for i = opt.steadycycle:cycle(end)-1,
        iscycle = cycle == i;
        
        b = last(iscycle) + 1;
        if (b > length(t))
            b = length(t);
        end;
        iscycle(b) = true;
        for j = 1:length(ang),
            axrot = S.Faxialtot(iscycle)*cos(ang(j)) - S.Flateraltot(iscycle)*sin(ang(j));
            
            Iax1(i,j) = trapz(t(iscycle),-axrot);
        end;
        a = find(sign(Iax1(i,1:end-1)) ~= sign(Iax1(i,2:end)));
        if (~isempty(a))
            b = a+1;
            ang0(i) = ang(a) + dang/(Iax1(i,b)-Iax1(i,a))*(0-Iax1(i,a));
        end;
    end;
    
    if (any(~isnan(ang0)))
        angrot = nanmedian(ang0);
        
        cosang = cos(angrot);
        sinang = sin(angrot);
        
        ax1 = S.Faxialt*cosang - S.Flateralt*sinang;
        lat1 = S.Faxialt*sinang + S.Flateralt*cosang;
        S.Faxialt = ax1;
        S.Flateralt = lat1;
        ax1 = S.Faxialn*cosang - S.Flateraln*sinang;
        lat1 = S.Faxialn*sinang + S.Flateraln*cosang;
        S.Faxialn = ax1;
        S.Flateraln = lat1;
        
        ax1 = S.Faxialtot*cosang - S.Flateraltot*sinang;
        lat1 = S.Faxialtot*sinang + S.Flateraltot*cosang;
        S.Faxialtot = ax1;
        S.Flateraltot = lat1;
        
        S.angrot = angrot;
    else
        warning('Couldn''t find a good angle)');
        S.angrot = NaN;
    end;
else
    S.angrot = 0;
end;

for i = 1:cycle(end),    
    iscycle = cycle == i;
   
    b = last(iscycle) + 1;
    if (b > length(t))
        b = length(t);
    end;
    iscycle(b) = true;
    %use last(iscycle)+1 to include the entire period.  Otherwise the sum
    %of all of the cycles wouldn't equal the integral over all time

    if (sum(iscycle) < 5)
        continue;
    end;
    
    for j = 1:size(S.Faxialt,1),
        %collect thrust due to normal and tangential fluid forces
        %note that S has reaction forces, so body force is the opposite
        tt = -S.Faxialt(j,iscycle,:);
        tt(tt < 0) = 0;
        tn = -S.Faxialn(j,iscycle,:);
        tn(tn < 0) = 0;
        
        %and integrate
        S.Ithrustt(j,i,1) = trapz(t(iscycle),tt(:,:,1));
        S.Ithrustn(j,i,1) = trapz(t(iscycle),tn(:,:,1));
        S.Ithrustt(j,i,2) = trapz(t(iscycle),tt(:,:,2));
        S.Ithrustn(j,i,2) = trapz(t(iscycle),tn(:,:,2));
        
        %if we were to do the same thing to get drag
        %  dn = S.Iaxialn(j,iscycle,:);
        %  dn(dn > 0) = 0;
        %  S.Idragn(j,i,1) = trapz(t(iscycle),dn(:,:,1));
        %then S.Ithrustn(j,i,1)+S.Idragn(j,i,1) would not equal
        %trapz(t(iscycle),S.Iaxialn(j,iscycle,1))
        %because of the way discrete integration works when we set certain
        %elements to zero.  So, instead we do the total integration and
        %subtract the thrust to get drag
        tott = trapz(t(iscycle),-S.Faxialt(j,iscycle,1));
        S.Idragt(j,i,1) = tott - S.Ithrustt(j,i,1);
        totn = trapz(t(iscycle),-S.Faxialn(j,iscycle,1));
        S.Idragn(j,i,1) = totn - S.Ithrustn(j,i,1);

        tott = trapz(t(iscycle),-S.Faxialt(j,iscycle,2));
        S.Idragt(j,i,2) = tott - S.Ithrustt(j,i,2);
        totn = trapz(t(iscycle),-S.Faxialn(j,iscycle,2));
        S.Idragn(j,i,2) = totn - S.Ithrustn(j,i,2);
    end;
    
    S.Iaxialtot(i) = trapz(t(iscycle),-S.Faxialtot(iscycle));
    S.Ilateraltot(i) = trapz(t(iscycle),-S.Flateraltot(iscycle));
end;
