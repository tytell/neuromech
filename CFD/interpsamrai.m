function varargout = interpsamrai(V, varargin)
% [x,y,V] = interpsamrai(V, ...)
%   or
% V = interpsamrai(V,x,y, ...)
%
% Options:
%   'numgrid' - 1 or 2 elements.  Specifies the number of points to span
%     the space represented in V.
%   'distgrid' - 1 or 2 elements.  Specifies the distance between points
%     that will span the space in V.
%   'smallest' - Generate a grid that spans the whole space with the
%     resolution of the finest grid in V.
%   'vars' - Specifies which variables to interpolate.
%
% Interpolates the multiscale grid in V to a single scale.

opt.numgrid = [];
opt.distgrid = [];
opt.smallest = false;
opt.vars = {};
opt.tol = 0.1;

isxy = false;
if ((nargin >= 3) && isnumeric(varargin{1}) && isnumeric(varargin{2}) && ...
        (ndims(varargin{1}) == ndims(varargin{2})) && ...
        all(size(varargin{1}) == size(varargin{2}))),
    [x,y] = deal(varargin{1:2});
    isxy = true;
    p = 3;
else
    p = 1;
end;
opt = parsevarargin(opt,varargin(p:end), p);

lo = cat(2,V.xlo);
hi = cat(2,V.xup);

if (~isxy),
    if (~isempty(opt.numgrid)),
        if (numel(opt.numgrid) == 1),
            opt.numgrid = opt.numgrid([1 1]);
        end;
        
        %generate the grid spanning the entire space
        [x,y] = meshgrid(linspace(min(lo(1,:)),max(hi(1,:)),opt.numgrid(1)), ...
            linspace(min(lo(2,:)),max(hi(2,:)),opt.numgrid(2)));
    elseif (~isempty(opt.distgrid)),
        if (numel(opt.distgrid) == 1),
            opt.distgrid = opt.distgrid([1 1]);
        end;
        
        xx = min(lo(1,:)):opt.distgrid(1):max(hi(1,:));
        xx = xx + (max(hi(1,:)) - xx(end))/2;
        yy = min(lo(2,:)):opt.distgrid(2):max(hi(2,:));
        yy = yy + (max(hi(2,:)) - yy(end))/2;
        
        [x,y] = meshgrid(xx,yy);
    elseif (opt.smallest),
        nc = cat(2,V.cols);
        nr = cat(2,V.rows);
        
        d = (hi(1:2,:) - lo(1:2,:)) ./ double([nc; nr]);
        d = min(d,[],2);
        
        xx = min(lo(1,:)):d(1):max(hi(1,:));
        xx = xx + (max(hi(1,:)) - xx(end))/2;
        yy = min(lo(2,:)):d(2):max(hi(2,:));
        yy = yy + (max(hi(2,:)) - yy(end))/2;
        
        [x,y] = meshgrid(xx,yy);
    end;
end;

if (isempty(opt.vars)),
    fieldnm = fieldnames(V);
else
    if (any(~ismember(opt.vars,fieldnames(V)))),
        error('Unknown var requested');
    end;
    fieldnm = opt.vars;
end;
nvar = length(fieldnm);
    

%set up our structure
vv = cell(1,nvar);
[vv{:}] = deal(NaN(size(x)));
V1 = cell2struct(vv,fieldnm,2);

%sort the patches by level (not usually necessary)
for i = 1:length(V),
    %interpolate only in the range of this patch
    inrng = (x >= lo(1,i)) & (x <= hi(1,i)) & ...
        (y >= lo(2,i)) & (y <= hi(2,i));
    
    if (any(inrng(:))),
        if (isfield(V,'dx')),
            dx = V(i).dx;
            dy = V(i).dy;
        else
            dx = (hi(1,i) - lo(1,i))/double(V(i).cols);
            dy = (hi(2,i) - lo(2,i))/double(V(i).rows);
        end;
        [x0,y0] = meshgrid(lo(1,i):dx:hi(1,i)+opt.tol*dx, ...
            lo(2,i):dy:hi(2,i)+opt.tol*dx);
        
        %interpolate each variable and put it into the structure,
        %overwriting the interpolated values from patches at lower levels
        for j = 1:nvar,
            if (all(size(V(i).(fieldnm{j})) == [V(i).rows V(i).cols])),
                v0 = V(i).(fieldnm{j});
                v0 = v0([1:end end],[1:end end]);
                v = interp2(x0,y0,v0, x(inrng),y(inrng));
                
                V1.(fieldnm{j})(inrng) = v;
            elseif (all(size(V(i).(fieldnm{j})) == [V(i).rows+1 V(i).cols+1])),
                v0 = V(i).(fieldnm{j});
                v = interp2(x0,y0,v0, x(inrng),y(inrng));
                
                if (any(flatten(isnan(v) & ~isnan(x(inrng))))),
                    beep;
                end;
                
                V1.(fieldnm{j})(inrng) = v;
            end;
        end;
    end;
end;

if ((nargout == 1) && isxy),
    varargout = {V1};
else
    varargout = {x,y,V1};
end;

