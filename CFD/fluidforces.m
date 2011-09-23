function varargout = fluidforces(varargin)
% F = fluidforces(t,freq, swimvecx,swimvecy, xm,ym,fxmtot,fymtot, xl,yl,fxltot,fyltot, ...
%    xr,yr,fxrtot,fyrtot)
%   or
% [Flateralt,Flateraln,Faxialt,Faxialn] = fluidforces(...)
%
% Options:
%  'nregions' - Number of regions along the body in which force will be
%     integrated.  The individual head and tail points are always included
%     in this number.
%  'debug" - Show debug plots
%  'ord' - Order of the tangent vector estimate

opt.ord = 4;
opt.nregions = 10;
opt.debug = false;
opt.dl = 4*pi/512;
opt.savevectors = false;
opt.saverawforce = false;
opt.method = 'wallstress';
opt.useallpoints = false;

if ((nargin == 0) || ((nargin > 0) && ischar(varargin{1}))),
    getvar('t','freq','xm','ym','fxmtot','fymtot','xl','yl','fxltot','fyltot', ...
        'xr','yr','fxrtot','fyrtot', ...
        'swimvecx','swimvecy');
    p = 1;
elseif ((nargin >= 6) && isstruct(varargin{5}) && isstruct(varargin{6})),
    [t,freq, swimvecx,swimvecy] = varargin{1:4};
    [IB,IBforce] = varargin{5:6};
    p = 6;
    use(IB);
    use(IBforce);
elseif ((nargin >= 16) && all(cellfun(@isnumeric,varargin(1:16))) && ...
        ((nargin == 16) || ischar(varargin{17})))
    [t,freq, swimvecx,swimvecy] = varargin{1:4};
    [xm,ym,fxmtot,fymtot, xl,yl,fxltot,fyltot, xr,yr,fxrtot,fyrtot] = varargin{5:16};
    p = 16;
elseif ((nargin >= 20) && all(cellfun(@isnumeric,varargin(1:20)))),
    [xm,ym,fxmtot,fymtot, xn,yn,fxntot,fyntot, xl,yl,fxltot,fyltot, xr,yr,fxrtot,fyrtot] = ...
        varargin{5:20};
    p = 20;
end;
    
opt = parsevarargin(opt,varargin(p+1:end), p);

npt = size(xl,1);
nfr = size(xl,2);

rgn = ones(2*npt-2,1);
rgn(1) = 1;
rgn(2:npt-1) = floor((1:npt-2)/(npt-1) * (opt.nregions-2)) + 2;
rgn(npt) = opt.nregions;
rgn(npt+1:end) = floor((npt-2:-1:1)/(npt-1) * (opt.nregions-2)) + 2;
side = repmat('L',[2*npt-2 1]);
side(1) = 'H';
side(npt) = 'T';
side(npt+1:end) = 'R';

[x,y, tanx,tany] = getoutlinetangent(xm,ym,xl,yl,xr,yr, 'ord',opt.ord);

fx = [fxmtot(1,:); fxltot(2:end-1,:); fxmtot(end,:); fxrtot(end-1:-1:2,:); NaN(1,nfr)];
fy = [fymtot(1,:); fyltot(2:end-1,:); fymtot(end,:); fyrtot(end-1:-1:2,:); NaN(1,nfr)];
if (opt.useallpoints)
    rgnmid = ones(2*npt-1,1);
    rgnmid(1) = 1;
    rgnmid(end) = opt.nregions;
    rgnmid(2:(2*npt-2)) = floor((1:(2*npt-3))/(2*npt-1) * (opt.nregions-2)) + 2;
    
    xmid = zeros(2*npt-1,nfr);
    xmid(1:2:end,:) = xm;
    xmid(2:2:end-1,:) = xn(1:end-1,:);
    ymid = zeros(2*npt-1,nfr);
    ymid(1:2:end,:) = ym;
    ymid(2:2:end-1,:) = yn(1:end-1,:);

    fxmid = zeros(2*npt-1,nfr);
    fxmid(1:2:end,:) = fxmtot;
    fxmid(2:2:end-1,:) = fxntot(1:end-1,:);
    fymid = zeros(2*npt-1,nfr);
    fymid(1:2:end,:) = fymtot;
    fymid(2:2:end-1,:) = fyntot(1:end-1,:);
    
    smid = [zeros(1,nfr); cumsum(sqrt(diff(xmid).^2 + diff(ymid).^2))];
    
    tanxmid = NaN(size(xmid));
    tanymid = NaN(size(ymid));
    tanxmid(2:end-1,:) = (xmid(3:end,:) - xmid(1:end-2,:)) ./ 2;
    tanymid(2:end-1,:) = (ymid(3:end,:) - ymid(1:end-2,:)) ./ 2;
    mag = sqrt(tanxmid.^2 + tanymid.^2);
    tanxmid = tanxmid ./ mag;
    tanymid = tanymid ./ mag;
    
    normxmid = tanymid;
    normymid = -tanxmid;
end;

swimvecx = repmat(swimvecx, [2*npt-1 1]);
swimvecy = repmat(swimvecy, [2*npt-1 1]);

normx = tany;
normy = -tanx;

s1 = [zeros(1,nfr); cumsum(sqrt(diff(x).^2 + diff(y).^2))];

switch lower(opt.method)
    case 'wallstress'
        dxds = deriv(s1(:,1),x);
        dyds = deriv(s1(:,1),y);
        dXds = sqrt(dxds.^2 + dyds.^2);
        
        %fluid stress based on Eq 8 and 9 in Williams et al 2009, Discrete and
        %Continuous Dynamical Systems series B 11(2): 519-540
        %these are NOT force densities -> they are actual forces
        fttot = -(fx.*tanx + fy.*tany) ./ dXds;
        fntot = -(fx.*normx + fy.*normy) ./ dXds;

        if (opt.useallpoints)
            dxmidds = deriv(smid(:,1),xmid);
            dymidds = deriv(smid(:,1),ymid);
            dXmidds = sqrt(dxmidds.^2 + dymidds.^2);
            
            ftmid = -(fxmid.*tanxmid + fymid.*tanymid) ./ dXmidds;
            fnmid = -(fxmid.*normxmid + fymid.*normymid) ./ dXmidds;
        end;
    case 'fluidstress'
        error('not implemented');
end;

faxialt = -fttot .* (tanx.*swimvecx + tany.*swimvecy);
faxialn = -fntot .* (normx.*swimvecx + normy.*swimvecy);
flateralt = -fttot .* (tanx.*swimvecy - tany.*swimvecx);
flateraln = -fntot .* (normx.*swimvecy - normy.*swimvecx);

if (opt.useallpoints)
    faxialtmid = - ftmid .* (tanxmid.*swimvecx + tanymid.*swimvecy);
    faxialnmid = - fnmid .* (normxmid.*swimvecx + normymid.*swimvecy);
    flateraltmid = - ftmid .* (tanxmid.*swimvecy - tanymid.*swimvecx);
    flateralnmid = - fnmid .* (normxmid.*swimvecy - normymid.*swimvecx);
end;

%integrate forces within each region on the body
%integrated, we get total forces -> N
if (opt.useallpoints)
    ns = 3;
else
    ns = 2;
end;
F.axialt = zeros(opt.nregions,nfr,ns);
F.axialn = zeros(opt.nregions,nfr,ns);
F.lateralt = zeros(opt.nregions,nfr,ns);
F.lateraln = zeros(opt.nregions,nfr,ns);
for j = 1:opt.nregions,
    isrgn = rgn == j;
    
    if (j == 1),
        %head is region one
        F.axialt(j,:,1) = faxialt(1,:);
    elseif (j == opt.nregions),
        %tail is the last region
        %same approximate "integration" procedure as the head
        F.axialt(j,:,1) = faxialt(npt,:);
    else
        for f = 1:nfr,
            %integrate each side for each frame separately
            good = isrgn & (side == 'L');
            F.axialt(j,f,1) = nansum(faxialt(good,f));
            F.axialn(j,f,1) = nansum(faxialn(good,f));
            F.lateralt(j,f,1) = nansum(flateralt(good,f));
            F.lateraln(j,f,1) = nansum(flateraln(good,f));
            
            good = isrgn & (side == 'R');
            F.axialt(j,f,2) = nansum(faxialt(good,f));
            F.axialn(j,f,2) = nansum(faxialn(good,f));
            F.lateralt(j,f,2) = nansum(flateralt(good,f));
            F.lateraln(j,f,2) = nansum(flateraln(good,f));
            
            if (opt.useallpoints)
                good = (rgnmid == j);
                F.axialt(j,f,3) = nansum(faxialtmid(good,f));
                F.axialn(j,f,3) = nansum(faxialnmid(good,f));
                F.lateralt(j,f,3) = nansum(flateraltmid(good,f));
                F.lateraln(j,f,3) = nansum(flateralnmid(good,f));
            end;
        end;
    end;
end;
F.axialtot = nansum(faxialt + faxialn,1);
F.lateraltot = nansum(flateralt + flateraln,1);
if (opt.useallpoints)
    F.axialtot = F.axialtot + nansum(faxialtmid + faxialnmid,1);
    F.lateraltot = F.lateraltot + nansum(flateraltmid + flateralnmid,1);
end;
    
%now we're going to calculate the total "thrust" and "drag" impulses per
%tail beat cycle.  We define "thrust" as axial force in the positive
%direction and "drag" as negative axial force
%integral gives impulse -> dynes sec 
cycle = floor(t*freq);
cycle = cycle - min(cycle) + 1;
ncycle = cycle(end);
F.thrustt = zeros(opt.nregions,ncycle,2);
F.thrustn = zeros(opt.nregions,ncycle,2);
F.dragt = zeros(opt.nregions,ncycle,2);
F.dragn = zeros(opt.nregions,ncycle,2);

for i = 1:cycle(end),    
    iscycle = cycle == i;
    if (sum(iscycle) <= 3)
        %skip anything with very few points - usually it's the last few
        %frames that might go into another cycle
        continue;
    end;
    
    for j = 1:opt.nregions,
        %collect thrust due to normal and tangential fluid forces
        tt = F.axialt(j,iscycle,:);
        tt(tt < 0) = 0;
        tn = F.axialn(j,iscycle,:);
        tn(tn < 0) = 0;
        
        %and integrate
        for k = 1:size(tt,3),
            F.thrustt(j,i,k) = trapz(t(iscycle),tt(:,:,k));
            F.thrustn(j,i,k) = trapz(t(iscycle),tn(:,:,k));
            
            %if we were to do the same thing to get drag
            %  dn = F.axialn(j,iscycle,:);
            %  dn(dn > 0) = 0;
            %  F.dragn(j,i,1) = trapz(t(iscycle),dn(:,:,1));
            %then F.thrustn(j,i,1)+F.dragn(j,i,1) would not equal
            %trapz(t(iscycle),F.axialn(j,iscycle,1))
            %because of the way discrete integration works when we set certain
            %elements to zero.  So, instead we do the total integration and
            %subtract the thrust to get drag
            tott = trapz(t(iscycle),F.axialt(j,iscycle,k));
            F.dragt(j,i,k) = tott - F.thrustt(j,i,k);
            
            totn = trapz(t(iscycle),F.axialn(j,iscycle,k));
            F.dragn(j,i,k) = totn - F.thrustn(j,i,k);
        end;
    end;
end;

if (opt.debug),
    %show some useful plots
    %total force is shown in a thick black line
    %left side is in green and right side is red
    %forces due to normal fluid force are solid lines, while those due to
    %tangential force are dotted
    if (getvar('t','freq')),
        phase = mod(t*freq,1);
        phase(diff(phase) < 0) = NaN;
        
        nper = floor(t(end)*freq);
        show = (t >= (nper-3)/freq) & (t <= nper/freq);
        
        Faxial = sum(F.axialt+F.axialn,3);
        Flateral = sum(F.lateralt+F.lateraln,3);
        
        a = [1 2 3 5 7 8 9 10];
        h = zeros(length(a),2);
        for i = 1:length(a),
            h(i,1) = subplot(length(a),2,2*i-1);
            ii = a(i);
            plot(phase(show), Faxial(ii,show), 'k-','LineWidth',2);
            addplot(phase(show),F.axialn(ii,show,1),'g-', phase(show),F.axialt(ii,show,1),'g:', ...
                phase(show),F.axialn(ii,show,2),'r-', phase(show),F.axialt(ii,show,2),'r:');
            ylabel(num2str(ii));
            
            h(i,2) = subplot(length(a),2,2*i);
            plot(phase(show), Flateral(ii,show), 'k-','LineWidth',2);
            addplot(phase(show),F.lateraln(ii,show,1),'g-', phase(show),F.lateralt(ii,show,1),'g:', ...
                phase(show),F.lateraln(ii,show,2),'r-', phase(show),F.lateralt(ii,show,2),'r:');    
        end;
        yl = [min(min(Faxial(:,show))) max(max(Faxial(:,show)))];
        set(h(:,1),'YLim',yl, 'XLim',[0 1]);
        yl = [min(min(Flateral(:,show))) max(max(Flateral(:,show)))];
        set(h(:,2),'YLim',yl, 'XLim',[0 1]);
        
        title(h(1,1),'Axial');
        title(h(1,2),'Lateral');
    end;
end;

if (opt.saverawforce),
    F.fttot = fttot;
    F.fntot = fntot;
end;

if (opt.savevectors),
    F.tanx = zeros(npt,nfr,2);
    F.tanx(:,:,1) = tanx(1:npt,:);
    F.tanx(:,:,2) = tanx(end:-1:npt,:);

    F.tany = zeros(npt,nfr,2);
    F.tany(:,:,1) = tany(1:npt,:);
    F.tany(:,:,2) = tany(end:-1:npt,:);
end;
F.rgn = rgn;

if (nargout == 0)
    putvar F;
elseif (nargout == 1)
    varargout = {F};
elseif (nargout == 4)
    varargout = {F.lateralt,F.lateraln,F.axialt,F.axialn};
elseif (nargout == 8),
    varargout = {F.lateralt,F.lateraln,F.axialt,F.axialn, ...
        flateralt,flateraln,faxialt,faxialn};
end;


