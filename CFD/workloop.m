function varargout = workloop(varargin)
% function [lennorm,worktot,workpos,workneg,workposact,worknegact,work] = workloop(t, x,y, fmus, act, opts...)
%   or       workloop(t, xl,yl,fmusl,actl, xr,yr,fmusr,actr, opts...)
%   or       workloop(t, xl,yl,ul,vl, fmusl,actl, xr,yr,ur,vr, fmusr,actr, opts...)

opt.ordlen = 15;
opt.plot = false;
opt.plotrange = [];
opt.plotpts = [0.25 0.5 0.75 0.95];
opt.musclescale = 38.6022;
opt.ps = 3;
opt.sidecol = 'br';
opt.margin = 0.1;
opt.gap = 0.03;
opt.per = 1;
opt.forcetaper = [];
opt.usevelocities = true;

[opt,args] = parsevarargin(opt, varargin, 1, 'allowno', 'leaveunknown');

if ((opt.usevelocities) && ((length(args) == 7) || (length(args) == 13)))
    [t, x,y,u,v, fmus,act] = deal(args{1:7});
    nsides = 1;
    if (length(args) == 13)
        [x2,y2,u2,v2, fmus2,act2] = deal(args{8:13});
        x = cat(3,x,x2);
        y = cat(3,y,y2);
        u = cat(3,u,u2);
        v = cat(3,v,v2);
        fmus = cat(3,fmus,fmus2);
        act = cat(3,act,act2);
        nsides = 2;
    end;
elseif ((length(args) == 5) || (length(args) == 9))
    [t, x,y, fmus,act] = deal(args{1:5});
    nsides = 1;
    if (length(args) == 9)
        [x2,y2, fmus2,act2] = deal(args{6:9});
        x = cat(3,x,x2);
        y = cat(3,y,y2);
        fmus = cat(3,fmus,fmus2);
        act = cat(3,act,act2);
        nsides = 2;
    end;
    opt.usevelocities = false;
end;

nfr = size(x,2);
npt = size(x,1);

if (opt.plot),
    if (isempty(opt.plotrange))
        opt.plotrange = [1 nfr];
        isplot = true(size(t));
    elseif (islogical(opt.plotrange) && (length(opt.plotrange) == length(t))),
        isplot = opt.plotrange;
        a = first(opt.plotrange);
        b = last(opt.plotrange);
        opt.plotrange = [a b];
    else
        isplot = false(size(t));
        isplot(opt.plotrange(1):opt.plotrange(2)) = true;
    end;
    
    nplot = length(opt.plotpts);

    height = (1-2*opt.margin - opt.gap*(nplot-1)) / nplot;
    
    h = -1*ones(nplot,3,nsides);
    
    width = 1-2*opt.margin - opt.gap - nsides*height;

    clf;
    
    width = width - opt.gap - height;
    hwork = axes('Position',...
        [1-opt.margin-height opt.margin ...
        height 1-2*opt.margin]);
    
    b = 1 - opt.margin - height;
    for i = 1:nplot,
        if (nsides == 1)
            height1 = height;
            b1 = b;
        else
            height1 = height/nsides;
            b1 = b + height1;
        end;
        h(i,1,1) = axes('Position',[opt.margin b1 width height1]);
        h(i,2,1) = axes('Position',[opt.margin b1 width height1], ...
            'Color','none');
        
        h(i,3,1) = axes('Position',[opt.margin+width+opt.gap b height height]);
        
        if (nsides == 2),
            b1 = b;
            h(i,1,2) = axes('Position',[opt.margin b1 width height1]);
            h(i,2,2) = axes('Position',[opt.margin b1 width height1], ...
                'Color','none');
            
            h(i,3,2) = axes('Position',...
                [opt.margin+width+opt.gap+height b height height]);
        end;
        
        b = b - height - opt.gap;
    end;
    
    lenrng = [Inf -Inf];
    Prng = [0 -Inf];
end;

nper = floor(range(t)/opt.per) + 1;
cycle = floor((t - t(1))/opt.per) + 1;

workpos = NaN(npt,nper,nsides);
workneg = NaN(npt,nper,nsides);
worktot = NaN(npt,nper,nsides);
workposact = NaN(npt,nper,nsides);
worknegact = NaN(npt,nper,nsides);

if (~isempty(opt.forcetaper)),
    fmus = fmus .* repmat(opt.forcetaper,[1 nfr nsides]);
end;
fmus = fmus/opt.musclescale * opt.ps;

if (opt.usevelocities)
    tanx = diff(x);
    tany = diff(y);
    mag = sqrt(tanx.^2 + tany.^2);
    mag(mag == 0) = 1;
    tanx = tanx ./ mag;
    tany = tany ./ mag;
    
    tanvel1 = u(1:end-1,:,:).*tanx + v(1:end-1,:,:).*tany;
    tanvel2 = u(2:end,:,:).*tanx + v(2:end,:,:).*tany;
    dlen = tanvel2 - tanvel1;
    dlen(npt,:,:) = 0;
    lennorm = [];
else
    %arclength
    s = [zeros(1,nfr,nsides); cumsum(sqrt(diff(x).^2 + diff(y).^2))];
    lenmn = nanmean(diff(s(:,1)));
    
    %length - calculate arc length between points separated by opt.ordlen
    len = zeros(size(s));
    orda = floor(opt.ordlen/2);
    ordb = ceil(opt.ordlen/2);
    len(orda+1:end-ordb,:,:) = (s(opt.ordlen+1:end,:,:) - s(1:end-opt.ordlen,:,:)) ./ ...
        opt.ordlen;
    len(len == 0) = lenmn;
    
    %normalize length according to the first frame
    %muscle force is in dynes
    lennorm = len./len(:,ones(1,nfr),ones(1,nsides)) - 1;
    
    %change in length in cm/s
    dlen = deriv(t,len,2);
end;

ispos = dlen < 0;

for i = 1:nper,
    iscycle = cycle == i;
    
    if (sum(iscycle) > 10)
        for side = 1:nsides,
            for pt = 1:npt,
                worktot(pt,i,side) = trapz(t(iscycle),-dlen(pt,iscycle,side) .* fmus(pt,iscycle,side));

                %NB: positive work is done when length is decreasing
                %integrate change in length * muscle force [cm/s] * [dynes]
                %over time [s] -> work in dynes * cm = ergs
                wp1 = -dlen(pt,iscycle,side) .* fmus(pt,iscycle,side);
                wp1(~ispos(pt,iscycle,side)) = 0;
                
                workpos(pt,i,side) = trapz(t(iscycle),wp1);
                workneg(pt,i,side) = worktot(pt,i,side) - workpos(pt,i,side);

                wa1 = -dlen(pt,iscycle,side) .* fmus(pt,iscycle,side);
                wa1(~act(pt,iscycle,side)) = 0;
                watot = trapz(t(iscycle), wa1);
                
                wpa1 = -dlen(pt,iscycle,side) .* fmus(pt,iscycle,side);
                wpa1(~act(pt,iscycle,side) | ~ispos(pt,iscycle,side)) = 0;

                workposact(pt,i,side) = trapz(t(iscycle), wpa1);
                worknegact(pt,i,side) = watot - workposact(pt,i,side);
            end;
        end;
    end;
end;
%instantaneous rate of working in dynes/cm * cm/s = ergs/s/cm
work = -dlen .* fmus;

if (opt.plot),
    pt = round(opt.plotpts*(npt-1)) + 1;

    for side = 1:nsides,
        col = opt.sidecol(side);

        for i = 1:nplot,
            plot(h(i,1,side), t(isplot), lennorm(pt(i),isplot), 'k:', 'LineWidth',2);
            
            xlim(h(i,1,side), t(opt.plotrange));
            if (side == 1),
                text('Parent',h(i,1,side), 'Units','normalized','Position',[0.01 1], ...
                    'String',sprintf('%c (%d%%)',char('A'+i-1),round(pt(i)/npt*100)), ...
                    'VerticalAlignment','top');
            end;
            if (nsides == 1),
                ylabel(h(i,1,side), 'length');
            else
                set(h(i,1:2,side), 'YTick',[], 'YColor','w');
            end;
            horizplot(h(i,1,side), 0,'Color',[0.7 0.7 0.7]);
            
            lenact = lennorm(pt(i),:);
            lenact(~act(pt(i),:)) = NaN;
            
            addplot(h(i,1,side), t(isplot), lenact(isplot), 'k-', 'LineWidth',3);
            
            fmusact = fmus(pt(i),:);
            fmusact(~act(pt(i),:)) = NaN;
            plot(h(i,2,side), t(isplot), fmus(pt(i),isplot), [col ':'], ...
                t(isplot),fmusact(isplot), [col '-'], 'LineWidth',2);
            xlim(h(i,2,side), t(opt.plotrange));
            set(h(i,2,side),'Color','none','YTickLabel',{});
            
            plot(h(i,3,side), lennorm(pt(i),isplot), fmus(pt(i),isplot), 'k-');
            addplot(h(i,3,side), lenact(isplot), fmus(pt(i),isplot), [col '-'], ...
                'LineWidth',3);
            
            set(h(i,3,side),'YAxisLocation','right');
            if (side == nsides),
                ylabel(h(i,3,side), 'force');
            end;
        end;
    end;
    
    lenrng1 = [min(flatten(lennorm(pt,isplot))) max(flatten(lennorm(pt,isplot)))];
    lenrng = [min(lenrng1(1),lenrng(1)) max(lenrng1(2),lenrng(2))];
    Prng1 = [0 max(flatten(fmus(pt,isplot)))];
    Prng = [0 max(Prng1(2),Prng(2))];
    
    pts = linspace(0,1,npt);
    cla(hwork,'reset');
    hold(hwork,'on');
    fill([0; nanmean(workposact(:,:,side),2); 0],pts([1 1:end end]), ...
        [0.7 0.7 0.7], 'EdgeColor','none', 'Parent',hwork);
    fill([0; nanmean(worknegact(:,:,side),2); 0],pts([1 1:end end]), ...
        [0.7 0.7 0.7], 'EdgeColor','none', 'Parent',hwork);
    addplot(hwork,nanmean(workpos(:,:,side),2),pts, [col '--'], ...
        nanmean(workneg(:,:,side),2),pts,[col '--']);
    addplot(hwork,nanmean(worktot(:,:,side),2),pts,[col '-'], 'LineWidth',3);
    addplot(hwork,[0 0],pts([1 end]),'k-');
    hold(hwork,'off');
end;

if (opt.plot),
    set(h(:,1,:),'YLim',lenrng);
    set(h(:,2,:),'YLim',Prng);
    linkaxes(flatten(h(:,1:2,:)),'x');
    
    set(h(:,3,:),'XLim',lenrng, 'YLim',Prng);
    linkaxes(flatten(h(:,3,:)),'xy');
    
    set(h(1:end-1,1:2,:),'XTick',[],'XColor','w');
    set(h(end,2,:), 'XTick',[], 'XColor','w');
    set(h(end,1:2,1:end-1), 'XTick',[],'XColor','w');
    
    if (nsides == 1),
        addscalebars(h(end,1,end),'x','xlen',0.5,'units','sec','Position',[1 -0.1]);
    else
        addscalebars(h(end,1,end),'xy','xlen',0.5, 'ylen',0.02, ...
            'units',{'sec',''}, 'Position',[-0.1 -0.1], 'Orientation','leftbottom');
    end;
    
    set(h(1:end-1,3,:), 'XTickLabel',{});
    if (nsides == 2),
        set(h(:,3,1), 'YTickLabel',{});
    end;
    
    xlabel(h(end,3), 'length');
    
    set(hwork,'YAxisLocation','right');
    horizplot(hwork, opt.plotpts, 'k-');
    axis(hwork,'ij','tight');
    xlabel(hwork,'work');
    ylabel(hwork,'position along body');
    varargout = {h,hwork};
else
    varargout = {lennorm,worktot,workpos,workneg,workposact,worknegact,work,dlen};
end;


