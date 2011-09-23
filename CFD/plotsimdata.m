function plotsimdata(file, varargin)
% function plotsimdata(file, varargin)
% Plots a summary of simulation data for the file specified.
% May save PDFs also.
% 'showfig' - indicates which figures to plot.  Can be a vector of values
%   from 1 to 8

opt.plotskip = 10;
opt.bodylen = 4*pi;
opt.savepdfs = false;
opt.showcycle = [];
opt.showfig = 1:8;

[opt] = parsevarargin(opt,varargin,1);

if (~isempty(file)),
    load(file);
    [pn,fn] = fileparts(file);
    load(fullfile(pn,[fn '_analysis.mat']));
end;

if (exist('good','var')),
    xm = xm(:,good);
    ym = ym(:,good);
    xn = xn(:,good);
    yn = yn(:,good);
    
    t = t(good);
    xl = xl(:,good);
    yl = yl(:,good);
    xr = xr(:,good);
    yr = yr(:,good);
    
    um = um(:,good);
    vm = vm(:,good);
    ul = ul(:,good);
    vl = vl(:,good);
    ur = ur(:,good);
    vr = vr(:,good);
    un = un(:,good);
    vn = vn(:,good);
    
    actl = actl(:,good);
    actr = actr(:,good);
    
    fmusPl = fmusPl(:,good);
    fmusPr = fmusPr(:,good);

    fxlmus = fxlmus(:,good);
    fylmus = fylmus(:,good);
    fxrmus = fxrmus(:,good);
    fyrmus = fyrmus(:,good);

    fxltot = fxltot(:,good);
    fyltot = fyltot(:,good);
    fxrtot = fxrtot(:,good);
    fyrtot = fyrtot(:,good);
    fxmtot = fxmtot(:,good);
    fymtot = fymtot(:,good);
    fxntot = fxntot(:,good);
    fyntot = fyntot(:,good);
end;

if (range(t(cyclenum == max(cyclenum))) < 0.99/freq),
    lastcycle = max(cyclenum)-1;
else
    lastcycle = max(cyclenum);
end;
issteady3 = issteady & (cyclenum >= lastcycle-3) & (cyclenum <= lastcycle);
steadypeak3 = true(size(indpeak));
good = isfinite(indpeak);
steadypeak3(good) = issteady3(indpeak(good));
steadypeak3 = any(steadypeak3);

if (isempty(opt.showcycle)),
    opt.showcycle = min(steadycycle + 2, lastcycle-1);
end;

steadycycles = steadycycle:size(indpeak,2);    

%% Table of values
if (ismember(1,opt.showfig)),
    figure(1);
    clf;
    
    ampmn = nanmean(amp(end,steadycycles));
    if (exist('flowspeed','var'))
        Umn = flowspeed;
        istether = true;
    else
        Umn = nanmean(comspeed(issteady));
        istether = false;
    end;
    Rebody = diground(Umn .* 4*pi / viscosity,100);
    Utail = nanmean(sqrt(um(end,issteady).^2 + vm(end,issteady).^2));
    widthtail = nanmean(width(end-2,:));
    Retail = diground(Utail .* widthtail / viscosity,1);
    
    actspeedmn = nanmean(actspeed(steadycycles));
    wavespeedmn = nanmean(wavespeed(:,steadycycles),2);
    wavelenmn = nanmean(flatten(wavelen(150:300,steadycycles)));
    
    St = diground(2*freq*ampmn / Umn,0.01);
    
    rownames = {'Simulation';'Frequency';'Amplitude';'Wave speed';...
        'Wave length';'Activ. speed';...
        'Wave/Activ';'Swim speed';'Re_{body}';'Re_{tail}';'St';'N steady cycles'};
    if (istether)
        rownames{8} = 'Flow speed';
    end;
    [~,fn] = fileparts(file);
    dat(:,1) = {fn; freq; ampmn / opt.bodylen; ...
        wavespeedmn/opt.bodylen; ...
        wavelenmn/opt.bodylen; ...
        actspeedmn/opt.bodylen; ...
        wavespeedmn/actspeedmn; ...
        Umn/opt.bodylen; ...
        Rebody;Retail; St; ...
        length(steadycycles)};
    dat(:,2) = {'';'Hz';'L';'L/s';'L';'L/s';'';'L/s';'';'';'';''};
    
    showtable(dat, 'rownames',rownames, 'align','lb', 'colgap',10);
    
    set(gcf,'Color','w');
end;

%% Plot simple kinematic parameters
% Body outlines, center of mass position, and swimming speed

if (ismember(2,opt.showfig)),
    figure(2);
    clf;
    
    %frames to show
    a = first(cyclenum == opt.showcycle);
    b = last(cyclenum == opt.showcycle)+1;
    ind = round(linspace(b,a,9));
    
    xx0 = comx(ind(5));
    yy0 = comy(ind(5));
    
    %estimate the head position if the speed was perfectly steady
    steadyposx = xm(1,ind(5)) + (t(ind)-t(ind(5)))*nanmean(comvelx(issteady));
    steadyposy = ym(1,ind(5)) + (t(ind)-t(ind(5)))*nanmean(comvely(issteady));
    
    subplot(5,1,1:2);
    
    svx = repmat(swimvecx,[npt 1]);
    svy = repmat(swimvecy,[npt 1]);
    tm = -(xm-xx0) .* svx - (ym-yy0) .* svy;
    nm = -(xm-xx0) .* svy + (ym-yy0) .* svx;
    tl = -(xl-xx0) .* svx - (yl-yy0) .* svy;
    nl = -(xl-xx0) .* svy + (yl-yy0) .* svx;
    tr = -(xr-xx0) .* svx - (yr-yy0) .* svy;
    nr = -(xr-xx0) .* svy + (yr-yy0) .* svx;
    steadypost = -(steadyposx-xx0) .* svx(1,ind) - (steadyposy-yy0) .* svy(1,ind);
    
    ot = [tl; tr(end:-1:1,:)];
    on = [nl; nr(end:-1:1,:)];
    fill(ot(:,ind)/opt.bodylen, (on(:,ind) + repmat(0:length(ind)-1,[2*npt 1]))/opt.bodylen, ...
        [0.7 0.7 0.7], 'EdgeColor','none');
    addquiverc(steadypost/opt.bodylen + 0.1, (0:length(ind)-1)/opt.bodylen, ...
        -0.1*ones(size(ind)), zeros(size(ind)), 'k', ...
        'AbsScale',1);
    
    for i = 1:length(ind),
        j = ind(i);
        
        isact = actl(:,j);
        istrans = [false; actl(1:end-1,j) ~= actl(2:end,j)];
        x1 = tl(:,j);
        y1 = nl(:,j);
        x1(istrans) = NaN;
        y1(istrans) = NaN;
        x1 = x1(isact | istrans);
        y1 = y1(isact | istrans);
        addplot(x1/opt.bodylen,(y1 + i-1)/opt.bodylen, 'b-', 'LineWidth',2);
        
        isact = actr(:,j);
        istrans = [false; actr(1:end-1,j) ~= actr(2:end,j)];
        x1 = tr(:,j);
        y1 = nr(:,j);
        x1(istrans) = NaN;
        y1(istrans) = NaN;
        x1 = x1(isact | istrans);
        y1 = y1(isact | istrans);
        addplot(x1/opt.bodylen,(y1 + i-1)/opt.bodylen, 'r-', 'LineWidth',2);
        
        for k = 1:size(indpeakcurve,2),
            if ((min(indpeakcurve(:,k)) <= j) && (max(indpeakcurve(:,k)) >= j)),
                c = find((indpeakcurve(1:end-1,k) <= j) & (indpeakcurve(2:end,k) > j));
                
                if (~isempty(c)),
                    p = round(c + (j-indpeakcurve(c,k))./(indpeakcurve(c+1,k)-indpeakcurve(c,k)));
                    
                    if (curve(p(1),j) > 0)
                        col = 'r';
                    else
                        col = 'b';
                    end;
                    addplot(tm(p,j)/opt.bodylen,(nm(p,j) + i-1)/opt.bodylen,[col '.']);
                end;
            end;
        end;
        
        for k = 1:size(indpeak,2),
            if ((min(indpeak(:,k)) <= j) && (max(indpeak(:,k)) >= j)),
                c = find((indpeak(1:end-1,k) <= j) & (indpeak(2:end,k) > j));
                
                if (~isempty(c)),
                    p = round(c + (j-indpeak(c,k))./(indpeak(c+1,k)-indpeak(c,k)));
                    
                    addplot(tm(p,j)/opt.bodylen,(nm(p,j) + i-1)/opt.bodylen,'k+');
                end;
            end;
        end;
    end;
    
    axis equal tight;
    axis off;
    
    addscalebars('x','xlen',0.1,'units','L');
    
    subplot(5,1,3:4);
    plot(s/opt.bodylen, ampcont(:,ind(1:end-1))/opt.bodylen, 'k-');
    axis equal tight;
    
    addscalebars('xy','position',[1.02 -0.05],'xlen',0.1,'ylen',0.1,'units',{'L','L'});
    axis off;
    
    subplot(5,1,5);
    plot(t,comspeed/opt.bodylen,'LineWidth',2);
    axis tight;
    
    xlabel('Time (s)');
    ylabel('Swimming speed');
    
    vertplot(first(t,issteady),'k:');
    horizplot(Umn/opt.bodylen,'k:');
    text(0,Umn/opt.bodylen, sprintf('%g L s-1',diground(Umn/opt.bodylen,0.01)), ...
        'VerticalAlignment','bottom');
    
    addscalebars('xy','units',{'sec','L s-1'})
    axis off;
    
    if (~isempty(file)),
        axes('Position',[0 0 1 1],'Units','normalized');
        text(0.5,0.99,file, 'HorizontalAlignment','center', 'VerticalAlignment','top');
        text(0.5,0.95,'Fig. 2. Kinematics', 'HorizontalAlignment','center',...
            'VerticalAlignment','top');
        axis off;
    end;
    
    set(gcf,'Color','w');
end;

%% Show curvature

if (ismember(3,opt.showfig)),
    figure(3);
    clf;
    
    k0 = [0.25 0.5 0.75 0.95];
    k = round(k0*npt);
    
    mplot(t(issteady3),curve(k,issteady3) .* opt.bodylen,...
        'k-|b--|g:|r-.','LineWidth',2);
    
    indpeak1 = indpeakcurve(k,steadypeak3);
    good = isfinite(indpeak1);
    tpk1 = NaN(size(indpeak1));
    tpk1(good) = t(indpeak1(good));
    cvpk1 = NaN(size(indpeak1));
    for i = 1:length(k),
        cvpk1(i,good(i,:)) = curve(k(i),indpeak1(i,good(i,:)));
    end;
    addmplot(tpk1',cvpk1' .* opt.bodylen,'ko|bs|gd|r^','MarkerSize',12, ...
        'MarkerFaceColor','w');
    
    grid;
    legend(num2str(k0'), ...
        'Location','NO','Orientation','horiz');
    legend boxoff;
    
    xlabel('Time (s)');
    ylabel('Curvature (L-1)');
    axis tight;
    xlim(t([first(issteady3) last(issteady3)]));
    
    if (~isempty(file)),
        axes('Position',[0 0 1 1],'Units','normalized');
        text(0.5,0.99,file, 'HorizontalAlignment','center', 'VerticalAlignment','top');
        text(0.5,0.95,'Fig. 3. Curvature', 'HorizontalAlignment','center',...
            'VerticalAlignment','top');
        axis off;
    end;
end;

iscurveleft = mode(sgnpeakcurve) < 0;
steadycycles = false(1,size(indpeak,2));
steadycycles(steadycycle:end) = true;

tpeak = NaN(size(indpeak));
good = isfinite(indpeak);
tpeak(good) = t(indpeak(good));

%% Show activation vs curvature
if (ismember(4,opt.showfig)),
    figure(4);
    clf;
    
    hwork = subplot(1,4,4);
    subplot(1,4,1:3);
    hold on;
    
    a = first(cyclenum == opt.showcycle);
    b = last(cyclenum == opt.showcycle)+1;
    
    isactcycle = any((indact <= b) & (indactoff >= a));
    
    phact = NaN(size(indact));
    good = isfinite(indact);
    phact(good) = t(indact(good)) * freq;
    phactoff = NaN(size(indactoff));
    good = isfinite(indactoff);
    phactoff(good) = t(indactoff(good)) * freq;
    
    phact0 = floor(min(phact));
    phact0 = matchsize(phact0,phact);
    phact = phact - phact0;
    phactoff = phactoff - phact0;
    
    ph1 = [phact(:,isactcycle); phactoff(end:-1:1,isactcycle)];
    ss = [s; s(end:-1:1)];
    isleft1 = isactleft(isactcycle);
    good = isfinite(ph1(:,1));
    lightblue = hsv2rgb([0.667 0.5 1]);
    h1 = fill(ph1(good,isleft1), ss(good)/opt.bodylen, lightblue,...
        'EdgeColor','b');
    good = isfinite(ph1(:,2));
    lightred = hsv2rgb([1 0.5 1]);
    h2 = fill(ph1(good,~isleft1), ss(good)/opt.bodylen, lightred,...
        'EdgeColor','r');
    
    phpeak = tpeak * freq;
    phpeakcurve = NaN(size(indpeakcurve));
    good = isfinite(indpeakcurve);
    phpeakcurve(good) = t(indpeakcurve(good)) * freq;
    
    phpeak = phpeak - phact0;
    phpeakcurve = phpeakcurve - phact0;
    
    h3 = addplot(phpeak(:,steadycycles & iscurveleft), s/opt.bodylen, 'b-');
    h4 = addplot(phpeak(:,steadycycles & ~iscurveleft), s/opt.bodylen, 'r-');
    h5 = addplot(phpeakcurve(:,steadycycles & iscurveleft), s/opt.bodylen, 'b:');
    h6 = addplot(phpeakcurve(:,steadycycles & ~iscurveleft), s/opt.bodylen, 'r:');
    
    hold off;
    
    ylabel('Position along body');
    xlabel('Time');
    axis tight;
    ylim([0 1]);
    legend([h1(1) h2(1) h3(1) h4(1) h5(1) h6(1)],'left act','right act',...
        'left exc', 'right exc', 'left curve', 'right curve',...
        'Location','NO','Orientation','horiz');
    legend boxoff;
    
    pos = get(gca,'Position');
    height = pos(4);
    pos = get(hwork,'Position');
    pos(4) = height;
    set(hwork,'Position',pos);
    
    if (~isempty(file)),
        axes('Position',[0 0 1 1],'Units','normalized');
        text(0.5,0.99,file, 'HorizontalAlignment','center', 'VerticalAlignment','top');
        text(0.5,0.95,'Fig. 4. Curvature and activation waves', 'HorizontalAlignment','center',...
            'VerticalAlignment','top');
        axis off;
    end;
end;

%% Show work loops

if (ismember(5,opt.showfig)),
    figure(5);
    clf;
    
    steady2 = last(issteady3);
    steady1 = first(issteady3);
    
    pts = [0.2 0.4 0.6 0.8 0.9];
    k = round(pts*(npt-1)) + 1;
    
    if (isforcetaper)
        forcetaper = width(:,1) ./ max(width(:,1));
    else
        forcetaper = ones(npt,1);
    end;
    
    [~,worktot,workpos,workneg,workposact,worknegact] = workloop(t,xr,yr,fmusPr,actr, xl,yl,fmusPl,actl, ...
        'plot',false,'per',1/freq,'forcetaper',forcetaper, 'ps',ps);
    hax = workloop(t,xr,yr,fmusPr,actr, 'plot', ...
        'plotrange',[steady1 steady2], 'sidecol','r', 'plotpts',pts, ...
        'per',1./freq,'forcetaper',forcetaper, 'ps',ps);
    
    for i = 1:length(pts),
        ispk = steadycycles & ~iscurveleft;
        vertplot(hax(i,1,1), tpeak(k(i),ispk), 'r-', 'LineWidth',2);
    end;

    if (~isempty(file)),
        axes('Position',[0 0 1 1],'Units','normalized');
        text(0.5,0.99,file, 'HorizontalAlignment','center', 'VerticalAlignment','top');
        text(0.5,0.95,'Fig. 5. Work loops', 'HorizontalAlignment','center',...
            'VerticalAlignment','top');
        axis off;
        set(gcf,'Color','w');
    end;
end;

if (ismember(4,opt.showfig)),
    snorm = linspace(0,1,npt)';
    iss = false(1,size(workpos,2));
    iss(steadycycle:end) = true;
    plot(hwork,workpos(:,iss,1),snorm, 'r:', ...
        workneg(:,iss,1),snorm, 'r--', ...
        workpos(:,iss,2),snorm, 'b:', ...
        workneg(:,iss,2),snorm, 'b--');
    hold(hwork,'on');
    fill([0; nanmean(workposact(:,:,1),2); 0],snorm([1 1:end end]), ...
        lightred, 'EdgeColor','none', 'Parent',hwork);
    fill([0; nanmean(worknegact(:,:,1),2); 0],snorm([1 1:end end]), ...
        lightred, 'EdgeColor','none', 'Parent',hwork);
    
    vertplot(hwork, 0, 'Color',[0.7 0.7 0.7]);
    
    addplot(hwork,worktot(:,iss,1),snorm, 'r-', ...
        worktot(:,iss,2),snorm, 'b-', ...
        'LineWidth',3);
    hold(hwork,'off');
    
    set(hwork,'YAxisLocation','right');
    xlabel(hwork,'work');
    ylabel(hwork,'Position along body');
    axis(hwork,'tight');
end;

if (ismember(6,opt.showfig)),
    figure(6);
    clf;
    
    %frames to show
    a = first(cyclenum == opt.showcycle);
    b = last(cyclenum == opt.showcycle)+1;
    ind = round(linspace(a,b,5));
    
    %rotate to point along the swimming direction
    svx = repmat(swimvecx,[npt 1]);
    svy = repmat(swimvecy,[npt 1]);
    tm = xm .* svx + ym .* svy;
    nm = xm .* svy - ym .* svx;
    tl = xl .* svx + yl .* svy;
    nl = xl .* svy - yl .* svx;
    tr = xr .* svx + yr .* svy;
    nr = xr .* svy - yr .* svx;
    
    %spacing between frames
    tspc = opt.bodylen * repmat(0:length(ind)-1, [2*npt 1]);
    nspc = zeros(2*npt, length(ind));
    
    h(1) = subplot(4,1,1);
    
    ot = [tl; tr(end:-1:1,:)];
    on = [nl; nr(end:-1:1,:)];
    fill(ot(:,ind)+tspc, on(:,ind) + nspc, ...
        [0.7 0.7 0.7], 'EdgeColor','none');
    
    pt = linspace(1,npt,size(Force.axialt,1)-1);
    pt = round((pt(1:end-1) + pt(2:end))/2);
    pt = [1 pt npt]';
    
    fax = sum(Force.axialt + Force.axialn, 3);
    flat = sum(Force.lateralt + Force.lateraln, 3);
    addquiverc(tm(pt,ind)+tspc(pt,:),nm(pt,ind) + nspc(pt,:), ...
        fax(:,ind),flat(:,ind),'k','s',0.8);
    axis equal tight;
    
    ylm = ylim;
    text(tm(pt(5),ind)+tspc(pt(5),:),ylm(2)*ones(size(ind)),num2str((1:length(ind))'), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom');
    axis off;
    
    h(2) = subplot(4,1,2);
    
    ph = mod(t*freq,1);
    ph([false diff(ph) < 0]) = NaN;
    
    plot(ph(issteady3),sum(Force.axialt(9,issteady3,:)+Force.axialn(9,issteady3,:),3),'k-',...
        'LineWidth',2);
    addplot(ph(issteady3),Force.axialt(9,issteady3,1),'b:', ph(issteady3),Force.axialn(9,issteady3,1),'b-', ...
        ph(issteady3),Force.axialt(9,issteady3,2),'r:', ph(issteady3),Force.axialn(9,issteady3,2),'r-');
    axis tight;
    vertplot(ph(ind),'k-');
    
    ylm = ylim;
    text(ph(ind),ylm(2)*ones(size(ind)),num2str((1:length(ind))'), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom');
    ylabel('Axial force');
    set(gca,'XTickLabel',{});
    
    h(3) = subplot(4,1,3);
    plot(ph(issteady3),sum(Force.lateralt(9,issteady3,:)+Force.lateraln(9,issteady3,:),3),'k-',...
        'LineWidth',2);
    addplot(ph(issteady3),Force.lateralt(9,issteady3,1),'b:', ph(issteady3),Force.lateraln(9,issteady3,1),'b-', ...
        ph(issteady3),Force.lateralt(9,issteady3,2),'r:', ph(issteady3),Force.lateraln(9,issteady3,2),'r-');    
    axis tight;
    vertplot(ph(ind),'k-');
    ylabel('Lateral force');
    xlabel('Phase');
    
    legend('total','L tan','L norm','R tan','R norm', ...
        'Location','EO','Orientation','vert');
    legend boxoff;

    h(4) = subplot(4,1,4);
    fthrust = sum(Force.thrustt(:,steadycycle:end,:) + Force.thrustn(:,steadycycle:end,:), 3);
    fdrag = sum(Force.dragt(:,steadycycle:end,:) + Force.dragn(:,steadycycle:end,:), 3);

    stairs([0; s(pt); s(end)]/opt.bodylen, fthrust([1 1:end end],:), 'g-');
    hold on;
    stairs([0; s(pt); s(end)]/opt.bodylen, fdrag([1 1:end end],:), 'r-');
    
    if (range(t(cyclenum == max(cyclenum))) < 0.99/freq),
        lastbeat = cyclenum == max(cyclenum)-1;
    else
        lastbeat = cyclenum == max(cyclenum);
    end;
    totalax = trapz(t(lastbeat),sum(Force.axialt(:,lastbeat,:) + Force.axialn(:,lastbeat,:), 3), 2);
    totalax = mean(totalax, 1);
    
    horizplot(totalax,'k-');
    text(1,totalax, sprintf('%.2f',totalax),'VerticalAlignment','middle');
    hold off;
    
    xlabel('Position along body');
    ylabel('Total impulse');
    axis tight;
    
    pos = get(h(3),'Position');
    w1 = pos(3);
    pos = get(h(2),'Position');
    h1 = pos(4);
    for i = 1:4,
        pos = get(h(i),'Position');
        pos(3) = w1;
        set(h(i),'Position',pos);
    end;
    pos = get(h(3),'Position');
    pos(4) = h1;
    set(h(3),'Position',pos);

    set(gcf,'Color','w');
    
    if (~isempty(file)),
        axes('Position',[0 0 1 1],'Units','normalized');
        text(0.5,0.99,file, 'HorizontalAlignment','center', 'VerticalAlignment','top');
        text(0.5,0.95,'Fig. 6. Wall stress force calculations', 'HorizontalAlignment','center',...
            'VerticalAlignment','top');
        axis off;
    end;
end;

IB = makeib(xl,yl,xr,yr,xm,ym,xn,yn,ul,vl,ur,vr,um,vm,un,vn);

if (ismember(7,opt.showfig)),
    figure(7);
    clf;
    
    Stress = fluidstressnearboundary(IB,{},swimvecx,swimvecy, 'fluidvals',Stress);
    
    if (~isfield(Stress,'Ithrustt'))
        Stress = stress2thrustdrag(t,cyclenum,Stress);
    end;
    
    %frames to show
    a = first(cyclenum == opt.showcycle);
    b = last(cyclenum == opt.showcycle)+1;
    ind = round(linspace(a,b,5));
    
    %rotate to point along the swimming direction
    svx = repmat(swimvecx,[npt 1]);
    svy = repmat(swimvecy,[npt 1]);
    tm = xm .* svx + ym .* svy;
    nm = xm .* svy - ym .* svx;
    tl = xl .* svx + yl .* svy;
    nl = xl .* svy - yl .* svx;
    tr = xr .* svx + yr .* svy;
    nr = xr .* svy - yr .* svx;
    
    %spacing between frames
    tspc = opt.bodylen * repmat(0:length(ind)-1, [2*npt 1]);
    nspc = zeros(2*npt, length(ind));
    
    h(1) = subplot(4,1,1);
    
    ot = [tl; tr(end:-1:1,:)];
    on = [nl; nr(end:-1:1,:)];
    fill(ot(:,ind)+tspc, on(:,ind) + nspc, ...
        [0.7 0.7 0.7], 'EdgeColor','none');
    
    edges = linspace(1,npt,size(Stress.Faxialt,1)+1);
    pt = round((edges(1:end-1) + edges(2:end))/2);
    
    fax = -sum(Stress.Faxialt + Stress.Faxialn, 3);
    flat = -sum(Stress.Flateralt + Stress.Flateraln, 3);
    addquiverc(tm(pt,ind)+tspc(pt,:),nm(pt,ind) + nspc(pt,:), ...
        fax(:,ind),flat(:,ind),'k','s',0.8);
    axis equal tight;
    
    yl = ylim;
    text(tm(pt(5),ind)+tspc(pt(5),:),yl(2)*ones(size(ind)),num2str((1:length(ind))'), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom');
    axis off;
    
    h(2) = subplot(4,1,2);
    
    ph = mod(t*freq,1);
    ph([false diff(ph) < 0]) = NaN;
    
    plot(ph(issteady3),-sum(Stress.Faxialt(end,issteady3,:)+Stress.Faxialn(end,issteady3,:),3),'k-',...
        'LineWidth',2);
    addplot(ph(issteady3),-Stress.Faxialt(end,issteady3,1),'b:', ph(issteady3),-Stress.Faxialn(end,issteady3,1),'b-', ...
        ph(issteady3),-Stress.Faxialt(end,issteady3,2),'r:', ph(issteady3),-Stress.Faxialn(end,issteady3,2),'r-');
    axis tight;
    vertplot(ph(ind),'k-');
    
    yl = ylim;
    text(ph(ind),yl(2)*ones(size(ind)),num2str((1:length(ind))'), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom');
    ylabel('Axial force');
    set(gca,'XTickLabel',{});
    
    h(3) = subplot(4,1,3);
    plot(ph(issteady3),-sum(Stress.Flateralt(end,issteady3,:)+Stress.Flateraln(end,issteady3,:),3),'k-',...
        'LineWidth',2);
    addplot(ph(issteady3),-Stress.Flateralt(end,issteady3,1),'b:', ph(issteady3),-Stress.Flateraln(end,issteady3,1),'b-', ...
        ph(issteady3),-Stress.Flateralt(end,issteady3,2),'r:', ph(issteady3),-Stress.Flateraln(end,issteady3,2),'r-');    
    axis tight;
    vertplot(ph(ind),'k-');
    ylabel('Lateral force');
    xlabel('Phase');
    
    legend('total','L tan','L norm','R tan','R norm', ...
        'Location','EO','Orientation','vert');
    legend boxoff;

    h(4) = subplot(4,1,4);
    Ithrust = sum(Stress.Ithrustt(:,steadycycle:lastcycle,:) + Stress.Ithrustn(:,steadycycle:lastcycle,:), 3);
    Idrag = sum(Stress.Idragt(:,steadycycle:lastcycle,:) + Stress.Idragn(:,steadycycle:lastcycle,:), 3);

    stairs(s(edges)/opt.bodylen, Ithrust([1:end end],:), 'g-');
    hold on;
    stairs(s(edges)/opt.bodylen, Idrag([1:end end],:), 'r-');
    
    if (range(t(cyclenum == max(cyclenum))) < 0.99/freq),
        lastbeat = cyclenum == max(cyclenum)-1;
    else
        lastbeat = cyclenum == max(cyclenum);
    end;
    totalax = mean(Stress.Iaxialtot(steadycycle:lastcycle));
    
    horizplot(totalax,'k-');
    text(1,totalax, sprintf('%.2f',totalax),'VerticalAlignment','middle');
    hold off;
    
    xlabel('Position along body');
    ylabel('Total impulse');
    axis tight;
    
    pos = get(h(3),'Position');
    w1 = pos(3);
    pos = get(h(2),'Position');
    h1 = pos(4);
    for i = 1:4,
        pos = get(h(i),'Position');
        pos(3) = w1;
        set(h(i),'Position',pos);
    end;
    pos = get(h(3),'Position');
    pos(4) = h1;
    set(h(3),'Position',pos);

    set(gcf,'Color','w');
    
    if (~isempty(file)),
        axes('Position',[0 0 1 1],'Units','normalized');
        text(0.5,0.99,file, 'HorizontalAlignment','center', 'VerticalAlignment','top');
        text(0.5,0.95,'Fig. 7. Fluid stress force calculations', 'HorizontalAlignment','center',...
            'VerticalAlignment','top');
        axis off;
    end;
end;
    
if (ismember(8,opt.showfig))
    Energy = energybalance(t,IB, ...
        makeibforce(fxltot,fyltot,fxrtot,fyrtot,fxmtot,fymtot,fxntot,fyntot, ...
        fxrmus,fyrmus,fxlmus,fylmus,fmusPr,fmusPl,actl,actr), swimvecx,swimvecy, {}, ...
        'fluidvals', Energy);

    wakekebycycle = Energy.totalkebycycle - Energy.bodykebycycle;
    ncycle = length(wakekebycycle);
    
    LHS = [Energy.axialkebycycle; Energy.lateralkebycycle; wakekebycycle; Energy.dissipbycycle];
    RHS = [Energy.springworkextbycycle; Energy.springworkintbycycle; Energy.muscleworkbycycle];
    
    figure(8);
    clf;

    cmap = jet(size(LHS,1)+size(RHS,1));
    barwidth = 0.3;
    
    LHSpatchx = zeros([4 ncycle size(LHS,1)]);
    LHSpatchy = zeros([4 ncycle size(LHS,1)]);
    LHScol = zeros(1,ncycle,size(LHS,1));
    barwidth1 = barwidth/size(LHS,1);
    for i = 1:size(LHS,1)
        bx1 = [(1:ncycle)-i*barwidth1; (1:ncycle)-i*barwidth1; ...
            (1:ncycle)-(i-1)*barwidth1; (1:ncycle)-(i-1)*barwidth1];
        by1 = [zeros(1,ncycle); LHS(i,:); LHS(i,:); zeros(1,ncycle)];
        
        LHSpatchx(:,:,i) = bx1;
        LHSpatchy(:,:,i) = by1;
        LHScol(1,:,i) = i;
    end;
    LHSsumx = [1:ncycle; 1:ncycle; (1:ncycle)-barwidth; (1:ncycle)-barwidth];
    LHSsumy = [zeros(1,ncycle); sum(LHS); sum(LHS); zeros(1,ncycle)];
    
    RHSpatchx = zeros([4 ncycle size(RHS,1)]);
    RHSpatchy = zeros([4 ncycle size(RHS,1)]);
    RHScol = zeros(1,ncycle,size(RHS,1));
    barwidth1 = barwidth/size(RHS,1);
    for i = 1:size(RHS,1),
        bx1 = [(1:ncycle)+(i-1)*barwidth1; (1:ncycle)+(i-1)*barwidth1; ...
            (1:ncycle)+i*barwidth1; (1:ncycle)+i*barwidth1];
        by1 = [zeros(1,ncycle); RHS(i,:); RHS(i,:); zeros(1,ncycle)];
        
        RHSpatchx(:,:,i) = bx1;
        RHSpatchy(:,:,i) = by1;
        RHScol(1,:,i) = i + size(LHS,1);
    end;
    RHSsumx = [1:ncycle; 1:ncycle; (1:ncycle)+barwidth; (1:ncycle)+barwidth];
    RHSsumy = [zeros(1,ncycle); sum(RHS); sum(RHS); zeros(1,ncycle)];
    
    clf;
    hold on;
    h3 = fill(LHSsumx,LHSsumy, 1, 'FaceColor','none','EdgeColor',cmap(1,:));
    h4 = fill(RHSsumx,RHSsumy, 1, 'FaceColor','none', 'EdgeColor',cmap(end,:));
    h1 = fill(flatten(LHSpatchx,2:3),flatten(LHSpatchy,2:3),flatten(LHScol,2:3), ...
        'EdgeColor','none');
    h2 = fill(flatten(RHSpatchx,2:3),flatten(RHSpatchy,2:3),flatten(RHScol,2:3), ...
        'EdgeColor','none');
    hold off;

    legend([h3(1); h1(1:ncycle:end); h4(1); h2(1:ncycle:end)], ...
        'LHS','    bodykeax','    bodykelat','    wakeke','    dissip',...
        'RHS','    extwork','    intwork','    muswork','Orientation','vertical', ...
        'Location','eastoutside');
    legend boxoff;
    
    xlabel('Tail beat cycle');
    ylabel('Change in energy');
    axis tight;
    xlim([0.5 ncycle+0.5]);
    
    set(gcf,'Color','w');
    
    if (~isempty(file)),
        axes('Position',[0 0 1 1],'Units','normalized');
        text(0.5,0.99,file, 'HorizontalAlignment','center', 'VerticalAlignment','top');
        text(0.5,0.95,'Fig. 8. Energy balance', 'HorizontalAlignment','center',...
            'VerticalAlignment','top');
        axis off;
    end;
end;

if (opt.savepdfs),
    [pn,fn] = fileparts(file);
    
    pdfoutname = fullfile(pn,sprintf('%s.pdf',fn));
    if (exist(pdfoutname,'file'))
        delete(pdfoutname);
    end;
    
    for f = opt.showfig,
        export_fig(pdfoutname, '-pdf','-append',figure(f));
    end;
end;



