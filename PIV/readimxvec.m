function varargout = readimxvec(file)
% function [x,y,u,v,w?,units] = readimxvec(file)

A = readimx(file);
if (~isfield(A,'Frames') || (length(A.Frames) ~= 1))
    error('Unrecognized IM7 file format');
end
A = A.Frames{1};

if ismember('W0',A.ComponentNames)
    nvec = 3;
else
    nvec = 2;
end

if any(~ismember({'U0','U1','U2','U3',...
        'V0','V1','V2','V3',...
        'ACTIVE_CHOICE','ENABLED'},A.ComponentNames))
    error('Unrecognized IM7 file format');
end

if (nvec == 3) && any(~ismember({'W0','W1','W2','W3'},A.ComponentNames))
    error('Unrecognized IM7 file format');
end

nms = genvarname(A.ComponentNames);
S = cell2struct(A.Components,nms,1);

u = S.U0.Planes{1};
v = S.V0.Planes{1};
if nvec == 3
    w = S.W0.Planes{1}
end

for i = 1:3
    txt = num2str(i);
    uname = ['U' txt];
    vname = ['V' txt];
    wname = ['W' txt];

    ischoice = S.ACTIVE_CHOICE.Planes{1} == i;
    
    u(ischoice) = S.(uname).Planes{1}(ischoice);
    v(ischoice) = S.(vname).Planes{1}(ischoice);
    if (nvec == 3)
        w(ischoice) = S.(wname).Planes{1}(ischoice);
    end
end

bad = (S.ACTIVE_CHOICE.Planes{1} == 4) | (S.ENABLED.Planes{1} == 0);
u(bad) = NaN;
v(bad) = NaN;
if (nvec == 3)
    w(bad) = NaN;
end

nx = size(u,1);
ny = size(u,2);

rx = 1:nx;
ry = 1:ny;

%set up the coordinate system
x = (rx-1)*A.Grids.X*A.Scales.X.Slope + A.Scales.X.Offset;
y = (ry-1)*A.Grids.Y*A.Scales.Y.Slope + A.Scales.Y.Offset;

%flip and transpose everything around so that we end up with normal
%Matlab order for coordinates
y = y(end:-1:1);

%make normal plaid matrices
[x,y] = meshgrid(x,y);

if isa(u,'single')
    u = double(u);
    v = double(v);
    if (nvec == 3)
        w = double(w);
    end
end

u = u'*A.Scales.I.Slope + A.Scales.I.Offset;
v = sign(A.Scales.Y.Slope)*v'*A.Scales.I.Slope + A.Scales.I.Offset;
u = flipud(u);
v = flipud(v);
if (nvec == 3),
    w = w'*A.Scales.I.Slope + A.Scales.I.Offset;
    w = flipud(w);
end;

units.x = A.Scales.X.Unit;
units.y = A.Scales.Y.Unit;
units.vel = A.Scales.I.Unit;

framedt = 1;
attrnames = cellfun(@(x) x.Name, A.Attributes, 'UniformOutput',false);
ind = find(strcmp(attrnames,'FrameDt'));
if (length(ind) == 1)
    txt = A.Attributes{ind}.Value;
    tok = regexp(txt,'([\d.]+)\s(\w*)','tokens','once');
    
    if (length(tok) == 2)
        framedt = str2double(tok{1});
        
        switch tok{2}
            case {'us','µs'}
                framedt = framedt/1e6;
            case 'ms'
                framedt = framedt/1000;
        end
    end
end
                
if (nvec == 3),
    varargout = {x,y,u,v,w};
else
    varargout = {x,y,u,v};
    if (nargout == 6),
        varargout{5} = [];
    end;
end;

if (nargout == length(varargout)+1),
    varargout{end+1} = units;
end;
