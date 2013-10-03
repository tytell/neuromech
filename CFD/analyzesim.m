function analyzesim(varargin)
% function analyzesim(filename,opts)
%   or     analyzesim()
% Runs the standard set of analysis procedures on a file, or in the base workspace.
% Calculates center of mass, swimming direction, kinematic parameters,
% muscle work, fluid stress and energy.  Saves to an output file
% ('outfile',filename).

opt.curvesmoothcutoff = 0.025;      % filter cutoff for curvature
% opt.nboundarypts = 30;
opt.steadythresh = 0.05;            % threshold velocity fluctuation.  Less than this is "steady"
opt.steadydX = 0.03;            % threshold velocity fluctuation.  Less than this is "steady"
opt.leftsidecurvature = -1;         % definition of "left" in terms of curvature
opt.leftsidepeak = 1;               % definition of "left" in terms of peak direction
% peak finding procedure
% Look for local max/min.  Find points > pkthresh*max.  Peak is center
% of all neighboring points > pkthresh*max.
opt.pkthresh = 0.9;                 
opt.outfile = '';                   % output file name
opt.samraibasedir = '';             % base directory for all the samrai data files
opt.rho = 1;                        % density in CGS units
opt.savepressure = false;
opt.savevorticity = false;

if ((nargin >= 1) && ischar(varargin{1}) && ...
        exist(varargin{1},'file')),
    datafilename = varargin{1};
    fprintf('**** File %s\n', varargin{1});
    filenameopt = {'-file',varargin{1}};
    p = 2;
else
    fprintf('Base workspace\n');
    datafilename = '';
    filenameopt = {};
    p = 1;
end;

opt = parsevarargin(opt, varargin(p:end), p);

% necessary variables
disp('  Loading variables...');
if (~getvar(filenameopt{:}, 'actl','actr','fmusPl','fmusPr',...
        'fxlmus','fxltot','fxmtot','fxntot','fxrmus','fxrtot',...
        'fylmus','fyltot','fymtot','fyntot','fyrmus','fyrtot',...        
        't','xl','yl','xm','ym','xn','yn','xr','yr',...
        'ul','vl','um','vm','un','vn','ur','vr'))
    error('Could not find necessary variables');
end;
if (getvar(filenameopt{:}, 'good')),
    actl = actl(:,good);
    actr = actr(:,good);
    fmusPl = fmusPl(:,good);
    fmusPr = fmusPr(:,good);
    fxlmus = fxlmus(:,good);
    fxrmus = fxrmus(:,good);
    fylmus = fylmus(:,good);
    fyrmus = fyrmus(:,good);
    fxmtot = fxmtot(:,good);
    fymtot = fymtot(:,good);
    fxntot = fxntot(:,good);
    fyntot = fyntot(:,good);
    fxltot = fxltot(:,good);
    fyltot = fyltot(:,good);
    fxrtot = fxrtot(:,good);
    fyrtot = fyrtot(:,good);
    t = t(good);
    xl = xl(:,good);
    yl = yl(:,good);
    xm = xm(:,good);
    ym = ym(:,good);
    xn = xn(:,good);
    yn = yn(:,good);
    xr = xr(:,good);
    yr = yr(:,good);

    ur = ur(:,good);
    vr = vr(:,good);
    um = um(:,good);
    vm = vm(:,good);
    un = un(:,good);
    vn = vn(:,good);
    ul = ul(:,good);
    vl = vl(:,good);
end;

%optional variables
%default values:
disp('  Getting the simulation parameters...');
freq = 1;
viscosity = 0.01;
sfo = 2.56e6;
sfo2 = '';
ps = 3;
gridres = 32;
isforcetaper = 'true';
if (~getvar(filenameopt{:}, 'freq','viscosity','sfo','sfo2','ps','gridres','isforcetaper',...
        '-keepundef')),
    prompt = {'freq','viscosity','sfo','sfo2','ps','gridres','isforcetaper'};
    def = cellfun(@num2str,{freq,viscosity,sfo,sfo2,ps,gridres,isforcetaper},'UniformOutput',false);
    if (isempty(sfo2)),
        def{4} = '0.2*2.56e6';
    end;
    if (isforcetaper),
        def{7} = 'true';
    else
        def{7} = 'false';
    end;
    vals = inputdlg(prompt,'Simulation constants',1,def);
    
    vals = cellfun(@eval,vals,'UniformOutput',false);
    
    [freq,viscosity,sfo,sfo2,ps,gridres,isforcetaper] = vals{:};
    
    putvar(filenameopt{:}, 'freq','viscosity','sfo','sfo2','ps','gridres','isforcetaper');
end;

nfr = size(xm,2);
npt = size(xm,1);

s0 = [zeros(1,nfr); cumsum(sqrt(diff(xm).^2 + diff(ym).^2))];
s = nanmean(s0,2);

istether = false;
disp('  Estimating com speed...');
if (all(abs(xm(1,:) - xm(1,1)) < 1e-3))
    disp('  -> appears to be a tethered head simulation');
    istether = true;
    
    comspeed = zeros(size(t));
    if (isempty(opt.samraibasedir))
        disp('  ... cannot get flow speed without samrai dirs');
        flowspeed = 0;
    else
        samraidirs = getdirectorynames(fullfile(opt.samraibasedir,'visit*'))';
        
        [x,y,V] = importsamrai(samraidirs{2},'vars',{'U_0','U_1'}, ...
            'interpolaten',[100 50]);
        flowspeed = mean(V.U_0(x < xm(1,1)-1));
    end;
end;

    
[width,area,sn, comx,comy,comvelx,comvely] = ...
    estcomspeed(t,xm,ym,xn,yn,xl,xr,yl,yr);

if (isforcetaper),
    forcetaper = width(:,1) ./ max(width(:,1));
else
    forcetaper = ones(npt,1);
end;

if (~istether)
    comspeed = sqrt(comvelx.^2 + comvely.^2);
    
    %unit vector pointing in the direction of the swimming speed, as a cubic
    %spline in 3 parts, which should get rid of any fluctuation at the level of
    %the tail beat frequency
    dt = t(2) - t(1);
    dur = round(2/dt);         %smooth over 2 sec
    comxsm = runavg(comx, dur);
    comysm = runavg(comy, dur);

    swimvecx = deriv(t,comxsm);
    swimvecy = deriv(t,comysm);
    mag = sqrt(swimvecx.^2 + swimvecy.^2);
    swimvecx = swimvecx ./ mag;
    swimvecy = swimvecy ./ mag;
else
    swimvecx = -1 * ones(size(t));
    swimvecy = zeros(size(t));
end;

%now average speeds over cycles to look for "steady"
%define cycles
cycle = t * freq;
cycle = cycle - min(floor(cycle));

cyclenum = floor(cycle)+1;
cyclenum2 = floor(cycle*2)+1;

disp('  Looking for steady speed...');
if (istether)
    ncycle = max(cyclenum);
    dt = t(2) - t(1);
    cycleindlen = round(1/freq/dt);
    
    k = 1:length(t)-cycleindlen;
    dX = NaN(size(t));
    dX(k+cycleindlen) = mean(sqrt((xm(:,k+cycleindlen) - xm(:,k)).^2 + ...
        (ym(:,k+cycleindlen) - ym(:,k)).^2));
    
    dXbycycle = accumarray(cyclenum(:), dX(:), [], @nanmean);
    steadycycle = length(dXbycycle) - 1;
    while ((steadycycle >= 1) && all(dXbycycle(steadycycle:end) < opt.steadydX))
        steadycycle = steadycycle-1;
    end;
else
    %average comspeed in each cycle
    comspeedbycycle = accumarray(cyclenum(:), comspeed(:), [], @nanmean);
    %and the amount it changes, relative to the current value
    dcomspeed = diff(comspeedbycycle) ./ comspeedbycycle(1:end-1);
    
    %look for those cycles that change less than opt.steadythresh
    steadycycle = length(comspeedbycycle)-1;
    while ((steadycycle >= 1) && all(dcomspeed(steadycycle:end) < opt.steadythresh)),
        steadycycle = steadycycle-1;
    end;
end;

if (steadycycle == 1)
    steadycycle = 2;
end
issteady = (cyclenum >= steadycycle);

%amplitude = displacement relative to the center of mass, perpendicular to the swimming
%vector
ampcont = -(xm - comx(ones(npt,1),:)) .* swimvecy(ones(npt,1),:) + ...
      (ym - comy(ones(npt,1),:)) .* swimvecx(ones(npt,1),:);
ampsteady = nanmean(range(ampcont(:,issteady)));
goodamp = issteady | (range(ampcont) < 2*ampsteady);
ampcont(:,~goodamp) = NaN;

%track curvature waves
%no smoothing in curvature, because we've already smoothed
disp('  Estimating curvature...');
curve = curvature(xm,ym, 'discrete','smooth',opt.curvesmoothcutoff);

%track zero crossings in curvature
disp('  Tracking curvature zero crossings...');
[indzero0,sgnzero0] = ...
    analyzeKinematics(s,t,xm,ym, 'curve',curve, 'nsmoothcurve',0, ...
                      'zeros', 'backwnd',0.5, 'returnpeaksign');
%track peaks in amplitude envelope
disp('  Tracking amplitude peaks...');
[indpeak0,sgnpeak0,~,~,~,~,~,wavespeed0,wavelen0,~,waven0] = ...
    analyzeKinematics(s,t,xm,ym, 'curve',ampcont, 'nsmoothcurve',10, ...
                      'peaks','backwnd',0.5, 'returnpeaksign');

%also track the activation wave
disp('  Finding the activation wave...');
[indact,indactoff,isactleft] = findactivationwave(actl,actr);

actlen = NaN(size(indact));
for j = 1:size(indact,2)-1,
    for i = 1:npt,
        if (isfinite(indact(i,j))),
            next = last(indact(:,j+1) <= indact(i,j));
            
            if (~isempty(next)),
                actlen(i,j) = 2*(s(i)-s(next));
            end;
        end;
    end;
end;

tact = NaN(size(indact));
good = isfinite(indact);
tact(good) = t(indact(good));

actspeed = NaN(1,size(tact,2));
for i = 1:length(actspeed),
    good = isfinite(tact(:,i));
    
    if (sum(good) > 100),
        p = polyfit(tact(good,i),s(good), 1);
    
        actspeed(i) = p(1);
    end;
end;

%early cycles are those with more than three activation starts on frame 1
early = sum(indact == 1) > 3;

firstright = first(~isactleft & ~early);
t0 = t(first(isfinite(indact(:,firstright))));

%match up the curvature to the activation
disp('  Aligning curvature and amplitude peaks to activation wave...');
cyclematchzero = zeros(1,size(indact,2));
cyclematchpeak = zeros(1,size(indact,2));
numincyclezero = zeros(1,size(indact,2));
numincyclepeak = zeros(1,size(indact,2));
for i = 1:size(indact,2),
    indact1 = first(indact(:,i),isfinite(indact(:,i)));
    indact2 = last(indactoff(:,i),isfinite(indactoff(:,i)));
    
    if (isactleft(i)),
        sgn1 = opt.leftsidecurvature;
    else
        sgn1 = -opt.leftsidecurvature;
    end;
    
    numincycle1 = sum((indzero0 >= indact1) & (indzero0 <= indact2));
    numincycle1(sgnzero0 ~= sgn1) = 0;
    
    [n1,a] = max(numincycle1);
    
    cyclematchzero(i) = a;
    numincyclezero(i) = n1;
    
    if (isactleft(i)),
        %NB: we want left side activation to correspond to peak right side
        %excursion, because that's appx when the left side muscle is
        %*shortest*
        sgn1 = -opt.leftsidepeak;
    else
        sgn1 = opt.leftsidepeak;
    end;
    
    numincycle1 = sum((indpeak0 >= indact1) & (indpeak0 <= indact2));
    numincycle1(sgnpeak0 ~= sgn1) = 0;
    
    [n1,a] = max(numincycle1);
    
    cyclematchpeak(i) = a;
    numincyclepeak(i) = n1;
end;

for i = 1:size(indzero0,2),
    k = find(cyclematchzero == i);
    if (length(k) > 1),
        [~,a] = max(numincyclezero(k));
        cyclematchzero(k) = 0;
        cyclematchzero(k(a)) = i;
    end;
end;

for i = 1:size(indpeak0,2),
    k = find(cyclematchpeak == i);
    if (length(k) > 1),
        [~,a] = max(numincyclepeak(k));
        cyclematchpeak(k) = 0;
        cyclematchpeak(k(a)) = i;
    end;
end;

%rearrange data to match the activation waves
indzero = NaN(size(indact));
good = cyclematchzero ~= 0;
indzero(:,good) = indzero0(:,cyclematchzero(good));
sgnzero = NaN(1,size(indact,2));
sgnzero(:,good) = sgnzero0(:,cyclematchzero(good));

indpeak = NaN(size(indact));
good = cyclematchpeak ~= 0;
indpeak(:,good) = indpeak0(:,cyclematchpeak(good));
sgnpeak = NaN(1,size(indact,2));
sgnpeak(good) = sgnpeak0(cyclematchpeak(good));

%find the peaks of curvature
disp('  Looking for curvature peaks...');
indpeakcurve = NaN(size(indzero));
sgnpeakcurve = NaN(size(indzero));
for j = 1:size(indzero,2)-1,
    for i = 1:size(indzero,1),
        if (~isnan(indzero(i,j)) && ~isnan(indzero(i,j+1)) && ...
                (indzero(i,j+1) > indzero(i,j))),
            k = indzero(i,j):indzero(i,j+1);
            curve1 = curve(i,k);
            mx = max(abs(curve1));
            
            ishigh = abs(curve1) >= opt.pkthresh*mx;
            sgn1 = sign(first(curve1,ishigh));
            if (any(sign(curve1(ishigh)) ~= sgn1)),
                ctr = NaN;
            else
                weight = abs(curve1(ishigh));
                if (length(weight) > 1),
                    weight = (weight - min(weight))/range(weight);
                else
                    weight = 1;
                end;
                ctr = sum(k(ishigh) .* weight) ./ ...
                      sum(weight);
            end;
            
            indpeakcurve(i,j) = round(ctr);
            sgnpeakcurve(i,j) = sgn1;
        end;
    end;
end;

good = cyclematchpeak ~= 0;

wavespeed = NaN(1,size(indzero,2));
wavespeed(good) = wavespeed0(cyclematchpeak(good));
waven = zeros(1,size(indzero,2));
waven(good) = waven0(cyclematchpeak(good));
wavelen = NaN(size(indzero));
wavelen(:,good) = wavelen0(:,cyclematchpeak(good));

amp = NaN(size(indact));
for c = 1:size(indact,2)-2,
    amp(:,c) = range(ampcont(:,(cyclenum2 == c) | (cyclenum2 == c+1)),2)/2;
end;

disp('  Calculating work loops...');
[lennorm,worktot,workpos,workneg,workposact,worknegact] = workloop(t, xl,yl,fmusPl,actl, ...
    xr,yr,fmusPr,actr, ...
    'per',1/freq, 'plot',false, 'forcetaper',forcetaper);

disp('  Calculating fluid forces...');
Force = fluidforces(t,freq, swimvecx,swimvecy, xm,ym, fxmtot,fymtot, ...
    xl,yl,fxltot,fyltot, xr,yr,fxrtot,fyrtot);

out.freq = freq;
out.viscosity = viscosity;
out.sfo = sfo;
out.sfo2 = sfo2;
out.ps = ps;
out.gridres = gridres;
out.isforcetaper = isforcetaper;
out.t = t;
out.s = s;
out.s0 = s0;
out.sn = sn;
out.nfr = nfr;
out.npt = npt;
out.comx = comx;
out.comy = comy;
out.comvelx = comvelx;
out.comvely = comvely;
out.swimvecx = swimvecx;
out.swimvecy = swimvecy;
out.comspeed = comspeed;
out.width = width;
out.area = area;
out.curve = curve;
out.ampcont = ampcont;
out.amp = amp;
out.wavespeed = wavespeed;
out.wavelen = wavelen;
out.indact = indact;
out.indactoff = indactoff;
out.isactleft = isactleft;
out.actspeed = actspeed;
out.actlen = actlen;
out.indzero = indzero;
out.sgnzero = sgnzero;
out.indpeak = indpeak;
out.sgnpeak = sgnpeak;
out.indpeakcurve = indpeakcurve;
out.sgnpeakcurve = sgnpeakcurve;
out.issteady = issteady;
out.steadycycle = steadycycle;
out.cyclenum = cyclenum;
out.lennorm = lennorm;
out.worktot = worktot;
out.workpos = workpos;
out.workneg = workneg;
out.workposact = workposact;
out.worknegact = worknegact;
out.Force = Force;
    
out.HGREV = savehgrev([], 'datafile',datafilename);


if (~isempty(opt.samraibasedir) && inputyn('Load fluid data? '))
    disp('  Calculating fluid quantities...');
    samraidirs = getdirectorynames(fullfile(opt.samraibasedir,'visit*'))';

    if (isempty(samraidirs))
        saveextra = {};
    else
        if (~isempty(opt.outfile))
            [pn,fn,ext] = fileparts(opt.outfile);
            
            cont = fullfile(pn,[fn 'Continue.mat']);
        else
            cont = '';
        end;
        
        if (~isempty(cont) && exist(cont,'file'))
            load(cont, 'Stress','samraidirs','fr','Energy');
            fr0 = fr;
            fprintf('Continuing from frame %d...\n', fr0);
        else
            fr0 = 1;
        end;
        
        boundx0 = [xm(1,:); xl(2:end-1,:); xm(end,:); xr(end-1:-1:2,:); xm(1,:)];
        boundy0 = [ym(1,:); yl(2:end-1,:); ym(end,:); yr(end-1:-1:2,:); ym(1,:)];
        boundu0 = [um(1,:); ul(2:end-1,:); um(end,:); ur(end-1:-1:2,:); um(1,:)];
        boundv0 = [vm(1,:); vl(2:end-1,:); vm(end,:); vr(end-1:-1:2,:); vm(1,:)];
        
        Energy.totalke = NaN(1,nfr);
        Energy.dissip = NaN(1,nfr);
        Energy.boundflux = NaN(4,nfr);
        
        opt1.rho = 1;               % in g/cm^3
        opt1.mu = viscosity;
        opt1.savepressure = opt.savepressure;
        opt1.savevorticity = opt.savevorticity;
        
        N = min(length(samraidirs)-1, nfr);
        for fr = fr0:N,
            fprintf('Importing %s (%d%%)...\n', samraidirs{fr+1}, round((fr+1)/nfr*100));
            V = importsamrai(samraidirs{fr+1},'vars',{'U_0','U_1','P'});
            V = getsamraipatchedges(V);
            V = velderivsamrai(V);
            
            [totalke1,dissip1,boundflux1] = energybalance1(V,opt1);
            Energy.totalke(fr) = totalke1;
            Energy.dissip(fr) = dissip1;
            Energy.boundflux(:,fr) = boundflux1;
            
            S1 = fluidstressnearboundary1(V, boundx0(:,fr),boundy0(:,fr),...
                boundu0(:,fr),boundv0(:,fr), opt1);
            
            Stress.tanx{fr} = S1.tanx;
            Stress.tany{fr} = S1.tany;
            Stress.s{fr} = S1.s;
            Stress.n{fr} = S1.n;
            Stress.nearboundx{fr} = S1.nearboundx;
            Stress.nearboundy{fr} = S1.nearboundy;
            Stress.normstress{fr} = S1.normstress;
            Stress.tanstress{fr} = S1.tanstress;
            if (opt.savepressure)
                Stress.nearboundp{fr} = S1.nearboundp;
            end;
            if (opt.savevorticity)
                Stress.nearboundut{fr} = S1.nearboundut;
                Stress.nearboundun{fr} = S1.nearboundun;
                Stress.vorticity{fr} = S1.vorticity;
            end;
            
            if (~isempty(cont) && (mod(fr,10) == 0))
                save(cont,'Stress','samraidirs','fr','Energy');
            end;
        end;
        
        out.Energy = Energy;
        out.Stress = Stress;
    end;
end;

if (istether)
    out.flowspeed = flowspeed;
end;

disp('  Saving data...');
if (~isempty(opt.outfile)),
    save(opt.outfile, '-struct','out');
else
    putvar(filenameopt{:}, '-fromstruct','out');
end;

    