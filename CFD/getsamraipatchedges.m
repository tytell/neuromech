function V = getsamraipatchedges(V)
% function V = getsamraipatchedges(V)
% Adds in the overlapping edges for each of the appropriate regions in V.
% Makes integration with trapzsamrai actually work properly.

opt.tol = 0.01;

npatch = length(V);

lo = cat(2,V.xlo);
hi = cat(2,V.xup);

for i = 1:npatch,
    V(i).dx = (V(i).xup(1) - V(i).xlo(1))/double(V(i).cols);
    V(i).dy = (V(i).xup(2) - V(i).xlo(2))/double(V(i).rows);
    V(i).px = lo(1,i):V(i).dx:hi(1,i)+opt.tol*V(i).dx;
    V(i).py = (lo(2,i):V(i).dy:hi(2,i)+opt.tol*V(i).dy)';
end;

vars = fieldnames(V)';
nvar = length(vars);
good = false(size(vars));
for i = 1:npatch,
    for v = 1:nvar,
        var1 = vars{v};
        if (size(V(i).(var1)) == [V(i).rows V(i).cols]),
            V(i).(var1) = V(i).(var1)([1:end end],[1:end end]);
            good(v) = true;
        end;
    end;
end;    
vars = vars(good);
nvar = length(vars);

for i = 1:npatch,
    ontop = find((abs(lo(2,:) - hi(2,i)) < opt.tol*V(i).dy) & ...
        (lo(1,:) < hi(1,i)) & (hi(1,:) > lo(1,i)));
    
    for j = ontop,
        %we have to do this stupid rounding because sometimes the px values
        %are not exactly equal, when they should be.  We know that they're
        %multiples of the same grid size, so we do the rounding.
        %Otherwise, sometimes the first or last interpolated value will be
        %NaN, because the first pxj value is 10^-18 or so greater than the
        %first pxi
        dx1 = min(V(i).dx, V(j).dx);
        pxi = round(V(i).px/dx1)*dx1;
        pxj = round(V(j).px/dx1)*dx1;
        
        inrng = (pxi >= lo(1,j)-opt.tol*dx1) & (pxi <= hi(1,j)-V(j).dx+opt.tol*dx1);
        inrng(end) = false;
        
        interpind = round((pxi(inrng) - pxj(1)) / V(j).dx)+1;
        
        for v = 1:nvar,
            var1 = vars{v};
            V(i).(var1)(end,inrng) = V(j).(var1)(1,interpind);
%             top1 = interp1(pxj(1:end-1),V(j).(var1)(1,1:end-1), pxi(inrng));
%             V(i).(var1)(end,inrng) = top1;
        end;
    end;
    
    onright = find((abs(lo(1,:) - hi(1,i)) < opt.tol*V(i).dx) & ...
        (lo(2,:) < hi(2,i)) & (hi(2,:) > lo(2,i)));
    for j = onright,
        %same rounding issue here
        dy1 = min(V(i).dy, V(j).dy);
        pyi = round(V(i).py/dy1)*dy1;
        pyj = round(V(j).py/dy1)*dy1;
        
        inrng = (pyi >= lo(2,j)-opt.tol*dy1) & (pyi <= hi(2,j)-V(j).dy+opt.tol*dy1);
        inrng(end) = false;
        
        interpind = round((pyi(inrng) - pyj(1)) / V(j).dy)+1;

        for v = 1:nvar,
            var1 = vars{v};
            V(i).(var1)(inrng,end) = V(j).(var1)(interpind,1);
            %right1 = interp1(pyj(1:end-1),V(j).(var1)(1:end-1,1), pyi(inrng));
            %V(i).(var1)(inrng,end) = right1;
        end;
    end;

    oncorner = find((abs(lo(1,:) - hi(1,i)) < opt.tol*V(i).dx) & ...
        (abs(lo(2,:) - hi(2,i)) < opt.tol*V(i).dy));
    for j = oncorner,
        for v = 1:nvar,
            V(i).(var1)(end,end) = V(j).(var1)(1,1);
        end;
    end;
end;
