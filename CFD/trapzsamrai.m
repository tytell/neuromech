function [val,integ] = trapzsamrai(V, fcn)
% [val,pergrid] = trapzsamrai(V, fcn)
% Integrates a function over all of the regions represented in V.  fcn is a
% function handle or a lambda function.  Make sure you call
% getsamraipatchedges.m before calling trapzsamrai.
%
% For example:
%  KE = trapzsamrai(V, @(x) (x.U_0.^2 + x.U_1.^2))
%   integrates velocity squared over the multiscale grid.

opt.tol = 0.01;

npatch = length(V);

dx = cat(2,V.dx);
dy = cat(2,V.dy);
lo = cat(2,V.xlo);
hi = cat(2,V.xup);

patchlev = cat(2,V.level_number);

integ = zeros(npatch,1);
for i = 1:npatch,
    v1 = double(feval(fcn,V(i)));
    
    [px1,py1] = meshgrid(V(i).px, V(i).py);
    integ1 = trapz2(px1,py1, v1);

    overlap = find((patchlev == V(i).level_number+1) & ...
        (lo(1,:) < hi(1,i)) & (hi(1,:) > lo(1,i)) & ...
        (lo(2,:) < hi(2,i)) & (hi(2,:) > lo(2,i)));
    
    sub = 0;
    for j = overlap,
        rgnx = (px1(1,:) > V(j).xlo(1)-opt.tol*dx(j)) & (px1(1,:) < V(j).xup(1)+opt.tol*dx(j));
        rgny = (py1(:,1) > V(j).xlo(2)-opt.tol*dy(j)) & (py1(:,1) < V(j).xup(2)+opt.tol*dy(j));
        
        if (any(rgnx) && any(rgny)),
            sub1 = trapz2(px1(rgny,rgnx),py1(rgny,rgnx), v1(rgny,rgnx));
            sub = sub + sub1;
        end;
    end;
    
    integ(i) = integ1 - sub;
end;    

val = sum(integ);
