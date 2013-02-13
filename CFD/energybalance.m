function E = energybalance(t, IB, F, swimvecx,swimvecy, samraidirs, varargin)
% function Energy = energybalance(t, IB, F, swimvecx,swimvecy, samraidirs,
%                            options...)
% Calculates the Energy structure, which contains the work done by the
% muscles, the kinetic energy of the fluid, etc.  The time intensive
% portion of the analysis is the importing of fluid data and integrating on
% the adaptive mesh, which is done in energybalance1.  If that's already
% been done, then this can take the fluid values as an option
% ('fluidvals',Energy), which makes the remaining calculations much faster.
%
% Change in fluid KE + fluid dissipation + KE flux across the boundaries
% should equal rate of muscle work + rate of spring ring
%
% See Fauci, L. J. and Peskin, C. S. (1988). A computational model of aquatic animal 
% locomotion. J. Comput. Phys. 77, 85-108.


opt.rho = 1;        % density of water in CGS units (g/ml)
opt.mu = 0.01;      % viscosity of water in CGS units
opt.dl = 4*pi/512;  % length scale factor
% opt.rng = [-0.5 0.5 -50 50];    
opt.debug = false;  % show debug plots
opt.forcetaper = true;      % force tapers toward tail?
opt.per = 1;                % cycle period
opt.checkworkintegral = false;      % check integrals as well as derivatives?
opt.fluidvals = struct;     % previously calculated fluid values (see energybalance1)
opt.bodyregions = 0:0.1:1;

opt = parsevarargin(opt,varargin, 4);

npt = size(IB.xl,1);
nfr = size(IB.xl,2);

if (~isempty(opt.fluidvals))
    E = opt.fluidvals;
end;

width = sqrt((IB.xl - IB.xr).^2 + (IB.yl - IB.yr).^2);

midx = (IB.xl + IB.xr)/2;
midy = (IB.yl + IB.yr)/2;

midu = (IB.ul + IB.ur)/2;
midv = (IB.vl + IB.vr)/2;

s = [zeros(1,nfr); cumsum(sqrt(diff(midx).^2 + diff(midy).^2))];

bke1 = 0.5*opt.rho*(midu.^2 + midv.^2);

bodyke = zeros(1,nfr);

for i = 1:nfr,
    bodyke(i) = trapz(s(:,i), width(:,i) .* bke1(:,i));
end;

%fxmtot etc are in dynes/height and um etc are in cm/s.  So workrate should be in
%ergs/s/cm
E.workrate = zeros(npt,nfr,4);
E.workrate(:,:,1) = IB.um .* F.fxmtot + IB.vm .* F.fymtot;
E.workrate(:,:,2) = IB.un .* F.fxntot + IB.vn .* F.fyntot;
E.workrate(:,:,3) = IB.ul .* F.fxltot + IB.vl .* F.fyltot;
E.workrate(:,:,4) = IB.ur .* F.fxrtot + IB.vr .* F.fyrtot;

E.workratemus = zeros(npt,nfr,4);
E.workratemus(:,:,3) = IB.ul .* F.fxlmus + IB.vl .* F.fylmus;
E.workratemus(:,:,4) = IB.ur .* F.fxrmus + IB.vr .* F.fyrmus;

E.workratesp = zeros(npt,nfr,4);
E.workratesp(:,:,1:2) = E.workrate(:,:,1:2);
E.workratesp(:,:,3) = IB.ul .* (F.fxltot - F.fxlmus) + IB.vl .* (F.fyltot - F.fylmus);
E.workratesp(:,:,4) = IB.ur .* (F.fxrtot - F.fxrmus) + IB.vr .* (F.fyrtot - F.fyrmus);

%workrate should be equal to workratemus + workratesp

E.totalwork = cumtrapz(t,sum(sum(E.workrate,1),3));

width = sqrt((IB.xl - IB.xr).^2 + (IB.yl - IB.yr).^2);

midx = (IB.xl + IB.xr)/2;
midy = (IB.yl + IB.yr)/2;

midu = (IB.ul + IB.ur)/2;
midv = (IB.vl + IB.vr)/2;

s = [zeros(1,nfr); cumsum(sqrt(diff(midx).^2 + diff(midy).^2))];

bke1 = 0.5*opt.rho*(midu.^2 + midv.^2);
keax1 = 0.5*opt.rho*(midu.*repmat(swimvecx,[npt 1]) + midv.*repmat(swimvecy,[npt 1])).^2;
kelat1 = 0.5*opt.rho*(midu.*repmat(swimvecy,[npt 1]) - midv.*repmat(swimvecx,[npt 1])).^2;

E.bodyke = zeros(1,nfr);
E.bodykeax = zeros(1,nfr);
E.bodykelat = zeros(1,nfr);
    
for i = 1:nfr,
    E.bodyke(i) = trapz(s(:,i), width(:,i) .* bke1(:,i));
    E.bodykeax(i) = trapz(s(:,i), width(:,i) .* keax1(:,i));
    E.bodykelat(i) = trapz(s(:,i), width(:,i) .* kelat1(:,i));
end;

if (isfield(E,'totalke'))
    fr0 = last(~isnan(E.totalke));
else
    E.totalke = NaN(1,nfr);
    E.dissip = NaN(1,nfr);
    E.boundflux = NaN(4,nfr);
end;

if (opt.debug),
    hax1 = gca;
    plot(hax1,t,sum(sum(E.workrate,1),3),'k-','LineWidth',2);
    
    dke = deriv(t,totalke);
    
    hke(1) = addplot(hax1,t,dke,'b-','LineWidth',2);
    hdis(1) = addplot(hax1,t,dissip,'r-','LineWidth',2);
    htot(1) = addplot(hax1,t,dke + dissip,'g-','LineWidth',2);
end;

for fr = fr0:length(samraidirs)-1,
    fprintf('Importing %s (%d%%)...\n', samraidirs{fr+1}, round((fr+1)/nfr*100));
    
    if (opt.checkworkintegral)
        V0 = importsamrai(samraidirs{fr+1},'vars',{'U_0','U_1','F0','F1'});
    else
        V0 = importsamrai(samraidirs{fr+1},'vars',{'U_0','U_1'});
    end;
    
    V0 = getsamraipatchedges(V0);
    V0 = velderivsamrai(V0);
    
    [totalke1,dissip1,boundflux1] = energybalance1(V0);
    
    E.totalke(fr) = totalke1;
    E.dissip(fr) = dissip1;
    E.boundflux(:,fr) = boundflux1;
    
    if (opt.debug),
        dke = deriv(t,totalke);
        
        set(hke,'YData',dke);
        set(hdis,'YData',dissip);
        set(htot,'YData',dke + dissip);
        
        drawnow;
        putvar totalke dissip boundflux;
        if (opt.checkworkintegral)
            putvar totalwork2 workrate2;
        end;
    end;
end;

cycles = floor(t / opt.per) + 1;
maxcycle = max(cycles);
if (range(t(cycles == maxcycle) >= 0.95*opt.per))
    ncycle = maxcycle;
else
    ncycle = maxcycle-1;
end;

E.totalkebycycle = zeros(1,ncycle);
E.dissipbycycle = zeros(1,ncycle);
E.totalworkbycycle = zeros(1,ncycle);
E.muscleworkbycycle = zeros(1,ncycle);
E.springworkintbycycle = zeros(1,ncycle);
E.springworkextbycycle = zeros(1,ncycle);

E.bodykebycycle = zeros(1,ncycle);
E.axialkebycycle = zeros(1,ncycle);
E.lateralkebycycle = zeros(1,ncycle);

if (~isempty(opt.bodyregions)),
    bodypts = round(opt.bodyregions * (npt-1)) + 1;
    nrgn = length(bodypts) - 1;
    E.muscleworkbyrgn = zeros(nrgn,ncycle);
    
    E.muscleactposworkbyrgn = zeros(nrgn,ncycle);
    E.muscleactnegworkbyrgn = zeros(nrgn,ncycle);
end;    

prevke = 0;
prevbke = 0;
prevaxke = 0;
prevlatke = 0;
for i = 1:ncycle,
    iscycle = cycles == i;

    ke1 = last(E.totalke,iscycle);
    E.totalkebycycle(i) = ke1 - prevke;
    prevke = ke1;
    
    E.dissipbycycle(i) = trapz(t(iscycle),E.dissip(iscycle));

    bke = last(E.bodyke,iscycle);
    E.bodykebycycle(i) = bke - prevbke;
    prevbke = bke;
    
    axke = last(E.bodykeax,iscycle);
    E.axialkebycycle(i) = axke - prevaxke;
    prevaxke = axke;
    
    latke = last(E.bodykelat,iscycle);
    E.lateralkebycycle(i) = latke - prevlatke;
    prevlatke = latke;
    
    E.totalworkbycycle(i) = trapz(t(iscycle),sum(sum(E.workrate(:,iscycle,:),1),3));
    E.muscleworkbycycle(i) = trapz(t(iscycle),sum(sum(E.workratemus(:,iscycle,:),1),3));
    if (~isempty(opt.bodyregions)),
        for j = 1:nrgn,
            isrgn = bodypts(j):bodypts(j+1)-1;
            
            workrgn = E.workratemus(isrgn,iscycle,:);
            E.muscleworkbyrgn(j,i) = trapz(t(iscycle),sum(sum(workrgn,1),3));
            
            workact = workrgn;
            actrgn = cat(3,false(length(isrgn),sum(iscycle),2),...
                F.actl(isrgn,iscycle), F.actr(isrgn,iscycle));

            workact(~actrgn) = 0;
            workposact = workact;
            workposact(workposact < 0) = 0;
            
            E.muscleactposworkbyrgn(j,i) = trapz(t(iscycle),sum(sum(workposact,1),3));
            
            workacttot = trapz(t(iscycle),sum(sum(workact,1),3));
            E.muscleactnegworkbyrgn(j,i) = workacttot - E.muscleactposworkbyrgn(j,i);
        end;
    end;
    
    E.springworkextbycycle(i) = trapz(t(iscycle),sum(sum(E.workratesp(:,iscycle,3:4),1),3));
    E.springworkintbycycle(i) = trapz(t(iscycle),sum(sum(E.workratesp(:,iscycle,1:2),1),3));
end;
    
if (opt.debug)
    wakekebycycle = E.totalkebycycle - E.bodykebycycle;
    
    LHS = [E.axialkebycycle; E.lateralkebycycle; wakekebycycle; E.dissipbycycle];
    RHS = [E.springworkextbycycle; E.springworkintbycycle; E.muscleworkbycycle];
    
    clf;
    hold on;
    bar((1:ncycle)-0.1, LHS', 0.2, 'stacked');
    bar((1:ncycle)+0.1, RHS(1,:)', 0.2, 'r', 'stacked');
    bar((1:ncycle)+0.1, [zeros(size(LHS)); RHS(2:3,:)]', 0.2, 'stacked');
    addplot((1:ncycle)+0.1, sum(RHS), 'ro','MarkerFaceColor','w');
    hold off;
end;




    