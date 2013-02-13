function S = fluidstressnearboundary(IB, samraidirs, swimvecx, swimvecy, varargin)
% function S = fluidstressnearboundary(IB, samraidirs, swimvecx, swimvecy, varargin)
% Estimates fluid stress near the immersed boundary IB on the basis of flow
% calculations in files in samraidirs.  The time consuming part of the
% calculation is loading the samrai data and calculating the normal and
% tangential stresses, which is done with a subroutine,
% fluidstressnearboundary1.  If these calculations are already completed,
% the stress can be passed as a structure option:
% fluidstressnearboundary(..., 'fluidvals',S);
%
% Units:
% S.normstress -> dynes/cm^2
% S.Faxialt, etc -> integral of stress ds -> [dynes/cm^2] * [cm] = dynes/cm

opt.debug = false;
opt.tol = 0.1;
opt.mu = 0.01;          % in Poise = g/cm/s
opt.rho = 1;            % g/cm
opt.nregions = 10;
opt.continuationfile = '';
opt.fluidvals = struct([]);
opt.normaldistance = 2;
opt.savepressure = false;
opt.savevorticity = false;
opt.frames = [];

opt = parsevarargin(opt,varargin, 5);

if (isempty(opt.frames))
    nfr = size(IB.xm,2);
    frames = 1:nfr;
else
    frames = opt.frames;
    nfr = length(frames);
end;

%boundary
boundx0 = [IB.xm(1,frames); IB.xl(2:end-1,frames); IB.xm(end,frames); IB.xr(end-1:-1:2,frames); IB.xm(1,frames)];
boundy0 = [IB.ym(1,frames); IB.yl(2:end-1,frames); IB.ym(end,frames); IB.yr(end-1:-1:2,frames); IB.ym(1,frames)];
boundu0 = [IB.um(1,frames); IB.ul(2:end-1,frames); IB.um(end,frames); IB.ur(end-1:-1:2,frames); IB.um(1,frames)];
boundv0 = [IB.vm(1,frames); IB.vl(2:end-1,frames); IB.vm(end,frames); IB.vr(end-1:-1:2,frames); IB.vm(1,frames)];

npt = size(IB.xm,1);

s0 = [zeros(1,nfr); cumsum(sqrt(diff(boundx0).^2 + diff(boundy0).^2))];
stail = s0(npt);

%set up the stress structure
S.tanx = cell(1,nfr);
S.tany = cell(1,nfr);
S.nearboundx = cell(1,nfr);
S.nearboundy = cell(1,nfr);
S.normstress = cell(1,nfr);
S.tanstress = cell(1,nfr);
if (opt.savepressure)
    S.nearboundp = cell(1,nfr);
end;
if (opt.savevorticity)
    S.nearboundut = cell(1,nfr);
    S.nearboundun = cell(1,nfr);
    S.vorticity = cell(1,nfr);
end;
S.s = cell(1,nfr);
S.n = cell(1,nfr);
S.axialstresst = cell(1,nfr);
S.lateralstresst = cell(1,nfr);
S.axialstressn = cell(1,nfr);
S.lateralstressn = cell(1,nfr);
S.Faxialt = zeros(opt.nregions,nfr,2);
S.Faxialn = zeros(opt.nregions,nfr,2);
S.Flateralt = zeros(opt.nregions,nfr,2);
S.Flateraln = zeros(opt.nregions,nfr,2);

if (opt.debug)
    getvar F;
    
    subplot(4,1,4);
    hnetfluidstress = line('XData',[],'YData',[],'Color','k');
    hnetwallstress = line('XData',[],'YData',[],'Color','r');
end;

if (~isempty(opt.fluidvals))
    S = joinstructfields(S,opt.fluidvals);
    i0 = last(~cellfun(@isempty,S.s));
    if (isempty(i0))
        i0 = 1;
    end;
elseif (~isempty(opt.continuationfile) && exist(opt.continuationfile,'file'))
    load(opt.continuationfile, 'S','fr','frames');
    i0 = find(frames == fr);
    fprintf('Continuing from frame %d...\n', i0);
else
    i0 = 1;
end;

if (~opt.debug && isempty(opt.fluidvals))
    timedWaitBar(0,'Calculating stress');
end;

if (~isempty(samraidirs))
    frames = frames((frames >= 1) & (frames <= length(samraidirs)-1));
    N = length(frames);
    %run the subroutine, if needed
    for i = i0:N,
        fr = frames(i);
        fprintf('Importing %s (%d%%)...\n', samraidirs{fr+1}, round((i+1)/N*100));
        V = importsamrai(samraidirs{fr+1},'vars',{'U_0','U_1','P'});
        V = getsamraipatchedges(V);
    
        S1 = fluidstressnearboundary1(V, boundx0(:,i),boundy0(:,i), ...
            boundu0(:,i),boundv0(:,i), opt);
        S.tanx{i} = S1.tanx;
        S.tany{i} = S1.tany;
        S.s{i} = S1.s;
        S.n{i} = S1.n;
        S.nearboundx{i} = S1.nearboundx;
        S.nearboundy{i} = S1.nearboundy;
        S.normstress{i} = S1.normstress;
        S.tanstress{i} = S1.tanstress;
        if (opt.savepressure)
            S.nearboundp{i} = S1.nearboundp;
        end;
        if (opt.savevorticity)
            S.nearboundut{i} = S1.nearboundut;
            S.nearboundun{i} = S1.nearboundun;
            S.vorticity{i} = S1.vorticity;
        end;
        
        if (nargout == 0),
            putvar S;
        end;
        
        if (opt.debug),
            ds = diff(s1);
            ds = ds([1 1:end]);
            
            subplot(4,1,1);
            plot(s0,F.fttot(:,i),'r-', s1,ds.*S.tanstress{i}(2,:),'k-');
            ylabel('tanforce');
            axis tight;
            
            subplot(4,1,2);
            plot(s0,F.fntot(:,i),'r-', s1,ds.*S.normstress{i}(2,:),'k-');
            ylabel('normforce');
            axis tight;
            
            subplot(4,1,3);
            cla reset;
            plot(boundx0(:,i),boundy0(:,i),'k-');
            
            k = 1:30:length(boundx);
            addquiverc(boundx(k),boundy(k), ...
                S.tanstress{i}(2,k).*tanx1(k) + S.normstress{i}(2,k).*normx1(k), ...
                S.tanstress{i}(2,k).*tany1(k) + S.normstress{i}(2,k).*normy1(k));
            axis equal tight;
            drawnow;
        else
            timedWaitBar(i/N);
        end;
        
        if (~isempty(opt.continuationfile) && (mod(i,10) == 0))
            save(opt.continuationfile,'S','samraidirs','fr','frames');
        end;
    end;
else
    N = nfr;
end;

if (numel(swimvecx) == 1)
    swimvecx = swimvecx*ones(1,N);
    swimvecy = swimvecy*ones(1,N);
end;


s0 = S.s{1}(1,:);
s0(s0 <= stail) = s0(s0 <= stail)/stail;
s0(s0 > stail) = (s0(s0 > stail) - stail) / (s0(end) - stail) + 1;

S.s0 = s0;
S.p = zeros(length(s0),nfr);

if (opt.savevorticity)
    S.nearboundvort = S.vorticity;
    S.vorticity = zeros(length(s0),nfr);
end;

%change coordinate systems from normal/tangential to axial/lateral and
%integrate along the regions of the body
for fr = 1:N,
    tanx1 = S.tanx{fr};
    tany1 = S.tany{fr};
    normx1 = tany1;
    normy1 = -tanx1;
    s1 = S.s{fr}(1,:);
    
    nind = opt.normaldistance;
    
    S.axialstresst{fr} = -S.tanstress{fr}(nind,:) .* (tanx1.*swimvecx(fr) + tany1.*swimvecy(fr));
    S.axialstressn{fr} = -S.normstress{fr}(nind,:) .* (normx1.*swimvecx(fr) + normy1.*swimvecy(fr));
    S.lateralstresst{fr} = -S.tanstress{fr}(nind,:) .* (tanx1.*swimvecy(fr) - tany1.*swimvecx(fr));
    S.lateralstressn{fr} = -S.normstress{fr}(nind,:) .* (normx1.*swimvecy(fr) - normy1.*swimvecx(fr));

    %find the tail tip
    dang = deriv(s1,tanx1).^2 + deriv(s1,tany1).^2;
    tail1 = first(s1 >= 0.95*stail);
    tail2 = last(s1 <= 1.05*stail);
    [~,pkind] = findpeaks2(dang(tail1:tail2),'max','sort','down');
    
    if (length(pkind) >= 2)
        %tail tip is the point midway between the two sharp corners
        stail1 = s1(round(tail1 + (pkind(1)+pkind(2))/2 - 1));
    else
        warning('Couldn''t find tail tip.  Using default...');
        stail1 = stail;
    end;
    
    rgnlen = s1(end)/((opt.nregions-1)*2);
    rgnlen2 = rgnlen/2;
    
    left1 = first(s1 > rgnlen2);
    left2 = last(s1 < stail1 - rgnlen2);
    
    isnose = (s1 <= rgnlen2) | (s1 >= s1(end)-rgnlen2);
    istail = (s1 >= stail1 - rgnlen2) & (s1 <= stail1 + rgnlen2);
    
    right1 = first(s1 > stail1 + rgnlen2);
    right2 = last(s1 < s1(end) - rgnlen2);
    
    rgn = zeros(size(s1));
    rgn(isnose) = 1;
    rgn(left1:left2) = floor((s1(left1:left2)-rgnlen2)/rgnlen)+2;
    rgn(istail) = opt.nregions;
    rgn(right1:right2) = floor((s1(end)-rgnlen2-s1(right1:right2))/rgnlen)+2;
    
    sidenum = zeros(size(s1));
    sidenum(isnose) = 1;
    sidenum(istail) = 1;
    sidenum(left1:left2) = 1;
    sidenum(right1:right2) = 2;
    
    %stresses have units of [g/cm/s^2] = [dynes/cm^2]
    %integrated over our boundary gives us [dynes/cm]
    %multiply by average width (0.29cm) to give an approximate 3D force, assuming
    %that we're kind of like a cylinder
    
    %integrate forces within each region on the body
    a = first(s1 >= s1(end)-rgnlen2);
    b = length(s1)-1;
    c = 1;
    d = last(s1 <= rgnlen2);
    
    noseind = [a:b c:d];
    xnose = S.nearboundx{fr}(1,noseind);
    ynose = S.nearboundy{fr}(1,noseind);
    snose = [0 cumsum(sqrt(diff(xnose).^2 + diff(ynose).^2))];

    S.Faxialt(1,fr,1) = trapz(snose,S.axialstresst{fr}(noseind));
    S.Faxialn(1,fr,1) = trapz(snose,S.axialstressn{fr}(noseind));
    S.Flateraln(1,fr,1) = trapz(snose,S.lateralstressn{fr}(noseind));
    S.Flateralt(1,fr,1) = trapz(snose,S.lateralstresst{fr}(noseind));
    
    for j = 2:opt.nregions,
        isrgn = rgn == j;
        
        for k = 1:2,
            %integrate each side for each frame separately
            a = first(isrgn & (sidenum == k));
            b = last(isrgn & (sidenum == k)) + 1;
            if (b > length(s1))
                %for the first region on the right side, which includes the
                %last point in the contour, the last point is equal to the
                %first point, so that the contour is closed.  So we don't have
                %to do anything special here to make sure that the integral
                %includes the entire wrapped contour
                b = length(s1);
            end;
            %use last + 1 so that the integral will include the entire contour.
            %Otherwise Faxialtot will not equal sum(sum(Faxialt+Faxialn,1),3)
            S.Faxialt(j,fr,k) = trapz(s1(a:b),S.axialstresst{fr}(a:b));
            S.Faxialn(j,fr,k) = trapz(s1(a:b),S.axialstressn{fr}(a:b));
            S.Flateraln(j,fr,k) = trapz(s1(a:b),S.lateralstressn{fr}(a:b));
            S.Flateralt(j,fr,k) = trapz(s1(a:b),S.lateralstresst{fr}(a:b));
        end;
    end;
    S.Faxialtot(fr) = trapz(s1,S.axialstresst{fr}+ ...
                            S.axialstressn{fr});
    
    if (opt.savepressure)
        ss = S.s{fr}(1,:);
        ss(ss <= stail1) = ss(ss <= stail1)/stail;
        ss(ss > stail1) = (ss(ss > stail1) - stail) / (S.s{fr}(1,end) - stail) + 1;
        
        S.p(:,fr) = interp1(ss',S.nearboundp{fr}(2,:)', s0');
    end;
    if (opt.savevorticity)
        ss = S.s{fr}(1,:);
        ss(ss <= stail1) = ss(ss <= stail1)/stail;
        ss(ss > stail1) = (ss(ss > stail1) - stail) / (S.s{fr}(1,end) - stail) + 1;
        
        S.vorticity(:,fr) = interp1(ss',S.nearboundvort{fr}(2,:)', s0');
    end;
    
    if (floor(10*(fr-1)/N) ~= floor(10*fr/N))
        fprintf('.');
    end;
end;
fprintf('\n');

if (isempty(opt.fluidvals))
  timedWaitBar(1);
end;

if (~isempty(opt.continuationfile))
    delete(opt.continuationfile);
end;

if (nargout == 0),
    putvar S;
end;
    
