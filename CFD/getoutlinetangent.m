function [xo,yo, tanx,tany] = getoutlinetangent(xm,ym, xl,yl, xr,yr, varargin)
% function [xo,yo, tanx,tany] = getoutlinetangent(xm,ym, xl,yl, xr,yr, varargin)
% Estimates the tangent vector to a set of IB points, accounting for the
% discontinuities near the head and tail properly.

opt.ord = 4;

opt = parsevarargin(opt, varargin, 6);

xo = [xm(1,:); xl(2:end-1,:); xm(end,:); xr(end-1:-1:2,:); xm(1,:)];
yo = [ym(1,:); yl(2:end-1,:); ym(end,:); yr(end-1:-1:2,:); ym(1,:)];

%calculate tangent vector based on distance between opt.ord points
orda = floor(opt.ord/2);
ordb = ceil(opt.ord/2);

xl([1 end],:) = xl([1 end],:);
yl([1 end],:) = yl([1 end],:);
xr([1 end],:) = xm([1 end],:);
yr([1 end],:) = ym([1 end],:);
xr = xr(end:-1:1,:);
yr = yr(end:-1:1,:);

tanxl = zeros(size(xl));
tanyl = zeros(size(yl));
tanxr = zeros(size(xr));
tanyr = zeros(size(yr));

%tangent vector at the head
tanxl(1,:) = (xl(ordb+1,:) - xr(end-orda,:)) / opt.ord;
tanyl(1,:) = (yl(ordb+1,:) - yr(end-orda,:)) / opt.ord;

%and the tail
tanxl(end,:) = (xr(ordb+1,:) - xl(end-orda,:)) / opt.ord;
tanyl(end,:) = (yr(ordb+1,:) - yl(end-orda,:)) / opt.ord;

%first points past the head and tail
tanxl(2,:) = (xl(3,:) - xm(1,:))/2;
tanyl(2,:) = (yl(3,:) - ym(1,:))/2;
tanxl(end-1,:) = (xm(end,:) - xl(end-2,:))/2;
tanyl(end-1,:) = (ym(end,:) - yl(end-2,:))/2;
tanxr(2,:) = (xr(3,:) - xm(end,:))/2;
tanyr(2,:) = (yr(3,:) - ym(end,:))/2;
tanxr(end-1,:) = (xm(1,:) - xr(end-2,:))/2;
tanyr(end-1,:) = (ym(1,:) - yr(end-2,:))/2;

%right sided difference for points close to the head
%on the left
tanxl(3:orda+1,:) = (xl(3+ordb:opt.ord+1,:) - xl(3:orda+1,:)) / ordb;
tanyl(3:orda+1,:) = (yl(3+ordb:opt.ord+1,:) - yl(3:orda+1,:)) / ordb;
%and close to the tail on the right (remember xr is now tail to head)
tanxr(3:orda+1,:) = (xr(3+ordb:opt.ord+1,:) - xr(3:orda+1,:)) / ordb;
tanyr(3:orda+1,:) = (yr(3+ordb:opt.ord+1,:) - yr(3:orda+1,:)) / ordb;

%left sided different for points close to the tail
%on the left
tanxl(end-ordb:end-2,:) = (xl(end-ordb:end-2,:) - xl(end-opt.ord:end-orda-2,:)) / orda;
tanyl(end-ordb:end-2,:) = (yl(end-ordb:end-2,:) - yl(end-opt.ord:end-orda-2,:)) / orda;
%and near the head on the right
tanxr(end-ordb:end-2,:) = (xr(end-ordb:end-2,:) - xr(end-opt.ord:end-orda-2,:)) / orda;
tanyr(end-ordb:end-2,:) = (yr(end-ordb:end-2,:) - yr(end-opt.ord:end-orda-2,:)) / orda;

%standard central different for middle points on both sides
tanxl(orda+2:end-ordb-1,:) = (xl(opt.ord+2:end-1,:) - xl(2:end-opt.ord-1,:)) ./ opt.ord;
tanyl(orda+2:end-ordb-1,:) = (yl(opt.ord+2:end-1,:) - yl(2:end-opt.ord-1,:)) ./ opt.ord;
tanxr(orda+2:end-ordb-1,:) = (xr(opt.ord+2:end-1,:) - xr(2:end-opt.ord-1,:)) ./ opt.ord;
tanyr(orda+2:end-ordb-1,:) = (yr(opt.ord+2:end-1,:) - yr(2:end-opt.ord-1,:)) ./ opt.ord;

tanx = cat(1,tanxl,tanxr(2:end-1,:),tanxl(1,:));
tany = cat(1,tanyl,tanyr(2:end-1,:),tanyl(1,:));

mag = sqrt(tanx.^2 + tany.^2);
tanx = tanx ./ mag;
tany = tany ./ mag;
