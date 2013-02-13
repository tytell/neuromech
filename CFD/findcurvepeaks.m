function [pospt,posfr] = findcurvepeaks(t,curve, per)
% function [pospt,posfr] = findcurvepeaks(t,curve, per)
% An attempt to put together a new algorithm for tracking peaks in
% curvature (different from the one in analyzekinematics), based on the
% fact that we know the period.  Not completed.

opt.nneighbor = 2;
opt.maxptjump = [-2 20];
opt.maxfrjump = [0 2];

npt = size(curve,1);

ispos = true(size(curve));
ispos(1:opt.nneighbor,:) = false;
ispos(end-opt.nneighbor+1:end,:) = false;
isneg = true(size(curve));
isneg(1:opt.nneighbor,:) = false;
isneg(end-opt.nneighbor+1:end,:) = false;

k = opt.nneighbor+1:npt-opt.nneighbor;
for off = 1:opt.nneighbor,
    ispos(k,:) = ispos(k,:) & (curve(k,:) > curve(k+off,:)) & ...
            (curve(k,:) > curve(k-off,:));
    isneg(k,:) = isneg(k,:) & (curve(k,:) < curve(k+off,:)) & ...
            (curve(k,:) < curve(k-off,:));
end;

[p,pfr] = find(ispos);
[n,nfr] = find(isneg);

plot(pfr, p, 'k.');

pospt = [];
posfr = [];
isconn = false(size(p));
while (any(~isconn)),
    next = first(~isconn);
    a = 1;

    addplot(pfr(next),p(next),'ro', 'MarkerFaceColor','r');

    pospt1 = [];
    posfr1 = [];
    
    hln = line('XData',posfr1, 'YData',pospt1, 'Color','r', 'Marker','o');
    hnext = line('XData',[], 'YData',[], 'Color','g', 'Marker','.', ...
                 'LineStyle','none');
    while (~isempty(next)),
        pospt1(a,1) = p(next);
        posfr1(a,1) = pfr(next);
        isconn(next) = true;
        
        set(hln,'XData',posfr1, 'YData',pospt1);
        
        dpt = p-pospt1(a);
        dfr = pfr-posfr1(a);
        
        a = a+1;
        
        dpt(isconn | (dpt < opt.maxptjump(1)) | (dpt > opt.maxptjump(2))) = NaN;
        dfr(isconn | (dfr < opt.maxfrjump(1)) | (dfr > opt.maxfrjump(2))) = NaN;
        
        if (any(~isnan(dpt) & ~isnan(dfr))),
            [d,next] = min(abs(dpt) + abs(dfr));
        
            set(hnext, 'XData',pfr(next), 'YData',p(next));
        
            drawnow;
        else
            next = [];
        end;
    end;
    
    addplot(posfr1,pospt1,'b.-');
    
    pospt = catuneven(2,pospt,pospt1);
    posfr = catuneven(2,posfr,posfr1);
end;

    
    
    

