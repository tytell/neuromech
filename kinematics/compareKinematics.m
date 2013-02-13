function compareKinematics(treat1files, treat2files)

m = 'os';
c = 'br';
c2 = 'cm';
off = 0.1;

treatments = ones(length(treat1files)+length(treat2files),1);
treatments(length(treat1files)+1:end) = 2;
files = {treat1files{:}, treat2files{:}};

indpeak1 = [];
for i = 1:length(files),
    load(files{i},'t','humms','hvmms');
    
    % rng(:,i) = [min(t) max(t)];
   
    clf;
    plot(t,humms, t,hvmms);
    xlabel('Time (s)');
    ylabel('Head velocity (mm/s)');
    legend('X','Y','Location','best');
    title(files{i});

    box = msgbox('Click time window to compare',['Set ' num2str(i)],'help','non-modal');

    yl = ylim;
    [tstart,q] = ginput(1);
    addplot([tstart tstart],yl,'k:');
    [tend,q] = ginput(1);
    addplot([tend tend],yl,'k:');
    
    if (ishandle(box)),
        delete(box);
    end;
    
    pause(0.5);

    rng(:,i) = [tstart; tend];
end;

colormap winter;
cmap = colormap;
t1colors = cmap(round(linspace(1,size(cmap,1), length(treat1files))),:);

colormap spring;
cmap = colormap;
t2colors = cmap(round(linspace(1,size(cmap,1), length(treat1files))),:);

colors = cat(1,t1colors,t2colors);
markers = repmat('o',size(treatments));
markers(treatments == 2) = 's';
linestyles = repmat('-',size(treatments));
linestyles(treatments == 2) = ':';

clf;
for i = 1:length(files),
    load(files{i},'t','humms','hvmms','hxmm','hymm','fishlenmm','indpeak',...
            'per','wavevel','wavelen','amp');
    linestyle1 = {linestyles(i), 'Color',colors(i,:)};
    switch treatments(i),
        case 1,
            dotstyle1 = {markers(i), 'MarkerFaceColor',colors(i,:), ...
                'MarkerEdgeColor','none','MarkerSize',6};
        case 2,
            dotstyle1 = {markers(i), 'Color',colors(i,:), 'MarkerSize',6};
    end;
    
    rngt = (t >= rng(1,i)) & (t <= rng(2,i));
    rngbt = (t(indpeak(end,:)) >= rng(1,i)) & (t(indpeak(end,:)) <= rng(2,i));
    peakind = indpeak(end,rngbt);
    
    nbeats(i) = sum(rngbt);
    
    %fit lines to the x and y head positions over time
    goodt = rngt & isfinite(hxmm(1,:)) & isfinite(hymm(1,:));
    px = polyfit(t(goodt),hxmm(1,goodt),1);
    py = polyfit(t(goodt),hymm(1,goodt),1);
    
    %get the general angle of travel
    ang(i) = atan2(py(1),px(1));
       
    %rotate head velocity vector
    par = humms*cos(ang(i)) + hvmms*sin(ang(i));
    parLs = par(peakind)/fishlenmm;
    perp = -humms*sin(ang(i)) + hvmms*cos(ang(i));
    
    swimvelLs(i) = nanmean(par(rngt))/fishlenmm;
    perpvelLs(i) = nanmean(perp(rngt))/fishlenmm;

    subplot(3,3,1);
    addplot(t(rngt),par(rngt)/fishlenmm,linestyle1{:});
    addplot(t(peakind),parLs,dotstyle1{:});
    ylabel('Swim vel (L/s)');
    xlabel('Time (s)');

    subplot(3,3,2);
    freq = 1./per(end,rngbt);
    addplot(parLs, freq, dotstyle1{:});
    ylabel('Tail freq (Hz)');

    subplot(3,3,3);
    stride = par(peakind)/fishlenmm.*per(end,rngbt);
    addplot(parLs,stride, dotstyle1{:});
    ylabel('Stride');
    
    subplot(3,3,4);
    addplot(parLs,wavevel(rngbt)/fishlenmm, dotstyle1{:});
    ylabel('Wave speed (L/s)');

    subplot(3,3,5);
    addplot(parLs,wavelen(end,rngbt)/fishlenmm,dotstyle1{:});
    ylabel('Wave len (L)');

    subplot(3,3,6);
    slip = par(peakind) ./ wavevel(rngbt);
    addplot(parLs,slip,dotstyle1{:});
    ylabel('Slip');
    
    subplot(3,3,7);
    addplot(parLs,amp(end,rngbt)/fishlenmm, dotstyle1{:});
    ylabel('Tail amp (L)');
    
    subplot(3,3,8);
    addplot(parLs,amp(1,rngbt)/fishlenmm, dotstyle1{:});
    ylabel('Head amp (L)');
    
    subplot(3,3,9);
    St = 2.*amp(end,rngbt) .* freq ./ par(indpeak(end,rngbt));
    addplot(parLs,St,dotstyle1{:});
    ylabel('St');
    xlabel('Swim vel (L/s)');
end;


