function h = h5orthoplot(h5file, varargin)
% function h = orthoplot(h5file)
%
% Plots volumetric data with three orthogonal slices like in the Zeiss
% LSM viewer program.  Allows the different slices to be moved with the
% mouse.  Takes x, y and z coordinates, plus the data itself in I.  The
% variable I must be three dimensional, but x,y, and z can be vectors or
% plaid matrices (like from meshgrid).
% Can also take a 3D vector field in vx,vy,vz, u,v,w and overlay that
% on top of the volumetric image.
%
% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2010, Eric Tytell

info = h5info(h5file,'/Image');
sz = info.Dataspace.Size;

attrnames = {info.Attributes.Name};
if ismember('Scale',attrnames)
    scale = h5readatt(h5file,'/Image','Scale');
else
    scale = 1;
end

x = ((1:sz(2)) - (sz(2)+1)/2) * scale;
y = ((1:sz(1)) - (sz(1)+1)/2) * scale;
z = ((1:sz(3)) - (sz(3)+1)/2) * scale;

% define the x and y positions of each orthogonal view
% since we're all in one axis, the xz and yz views have to be offset
% down and right, respectively.
xyx = x([1 end]);
xyy = y([1 end]);
xzx = x([1 end]);
xzy = [0 range(z)] + y(end);
yzx = [0 range(z)] + x(end);
yzy = y([1 end]);

% get the middle index from each coordinate array
% TODO: in principle, we shouldn't assume even spacing, as we do here
mid = ceil([length(x)/2 length(y)/2 length(z)/2]);

Ixy = h5read(h5file,'/Image',[1 1 mid(3)],[sz(1) sz(2) 1]);
Ixz = h5read(h5file,'/Image',[mid(2) 1 1],[1 sz(2) sz(3)]);
Iyz = h5read(h5file,'/Image',[1 mid(1) 1],[sz(1) 1 sz(3)]);

% plot the images
himxy = imagesc(xyx,xyy, Ixy);
hold on;
himxz = imagesc(xzx,xzy, squeeze(Ixz)');
himyz = imagesc(yzx,yzy, squeeze(Iyz));

% yes, the "tight" needs to be there twice.  Don't ask me why, but it does
% what I want it to that way
axis tight equal tight off;
colormap gray;

% draw the vectors
if (nargin == 10),
    % unplaid the x,y,z
    if ((ndims(vx) == ndims(u)) && all(size(vx) == size(u))),
        vx = vx(1,:,1);
        vy = vy(:,1,1);
        vz = squeeze(vz(1,1,:));
    end;
    % and find the middle index, taking into account possible NaNs
    [q,mx] = min(abs(vx-(x(end)+x(1))/2));
    [q,my] = min(abs(vy-(y(end)+y(1))/2));
    [q,mz] = min(abs(vz-(z(end)+z(1))/2));

    mid = [mx my mz];

    % offset for the xz and yz axes
    xzOffset = max(y) - min(z);
    yzOffset = max(x) - min(z);

    % calculate a 3D scaling value
    mag = sqrt(u.^2 + v.^2 + w.^2);
    if (numel(vx) > 1),
        d = nanmean(cat(1,diff(vx(:)),diff(vy(:)),diff(vz(:))));
    else
        d = min(size(I))/2;
    end;
    vscale = d/prctile(mag(:),95);

    % and draw the vectors
    hvecxy = quiverc(vx,vy, ...
                     u(:,:,mid(3)),v(:,:,mid(3)),'y',...
                     'AbsScale',vscale);
    hvecxz = quiverc(vx,vz+xzOffset, ...
                     squeeze(u(mid(2),:,:)),...
                     squeeze(w(mid(2),:,:)),'r','AbsScale',vscale);
    hvecyz = quiverc(vz+yzOffset,vy,...
                     squeeze(w(:,mid(1),:)),...
                     squeeze(v(:,mid(1),:)),'g','AbsScale',vscale);

    % set up the Vec structure
    data.Vec.x = vx;
    data.Vec.y = vy;
    data.Vec.z = vz;
    data.Vec.u = u;
    data.Vec.v = v;
    data.Vec.w = w;
    data.Vec.xzOffset = xzOffset;
    data.Vec.yzOffset = yzOffset;
    data.Vec.scale = vscale;
    data.Vec.h = {hvecxy,hvecxz,hvecyz};
end;

% box in the different views
plot(xyx([1 2 2 1 1]),xyy([1 1 2 2 1]),'y-',...
     xzx([1 2 2 1 1]),xzy([1 1 2 2 1]),'r-',...
     yzx([1 2 2 1 1]),yzy([1 1 2 2 1]),'g-',...
     'LineWidth',1.5);

% and draw the lines that indicate where the planes are
hline = repmat(-1,[3 2]);
hline(1,1) = plot([xyx(1) yzx(2)],mean(xyy)*[1 1],'r',...
                  'ButtonDownFcn',{@orthoLineDown,1});
hline(2,1) = plot(mean(xyx)*[1 1],[xyy(1) xzy(2)],'g',...
                  'ButtonDownFcn',{@orthoLineDown,2});
hline(3,1) = plot(xyx,mean(xzy)*[1 1],'y',...
                  'ButtonDownFcn',{@orthoLineDown,3});
hline(3,2) = plot(mean(yzx)*[1 1],yzy,'y',...
                  'ButtonDownFcn',{@orthoLineDown,3});

hold off;

% set up the guidata for the button handler routines
data.Figure = gcf;
set(data.Figure, 'DoubleBuffer', 'on', 'KeyPressFcn',@orthoKeyPress);

data.X = xyx;
data.Y = xyy;
data.Z = {xzy yzx};
data.indx = mid(1);
data.indy = mid(2);
data.indz = mid(3);
data.h5file = h5file;
data.sz = sz;

data.hImages = [himxy himxz himyz];
data.hLines = hline;
guidata(data.Figure, data);

% -------------------------------------------------------------
function data = orthoUpdate(data, indx,indy,indz)

if ~isempty(indy)
    data.indy = indy;
    
    pos = data.Y(1) + (indy - 1)/(data.sz(1)-1) * diff(data.Y);
    set(data.hLines(1,1),'YData',[pos pos]);
    
    Ixz = h5read(data.h5file,'/Image',[indy 1 1],[1 data.sz(2) data.sz(3)]);
    set(data.hImages(2), 'CData', squeeze(Ixz)');

    if (isfield(data,'Vec') && ...
            (pos >= data.Vec.y(1)) && (pos <= data.Vec.y(end))),
        vindy = (indy - 1)/(data.sz(1)-1) * (size(data.Vec.u,1)-1) + 1;

        delete(data.Vec.h{2});
        
        hold on;
        data.Vec.h{2} = quiverc(data.Vec.x,data.Vec.z+data.Vec.xzOffset, ...
            data.Vec.u(vindy,:,:),data.Vec.w(vindy,:,:), ...
            'AbsScale',data.Vec.scale,'r');
        hold off;
    end;
end;

if ~isempty(indx)
    data.indx = indx;
    pos = data.X(1) + (indx - 1)/(data.sz(2)-1) * diff(data.X);
    set(data.hLines(2,1),'XData',[pos pos]);
    
    Iyz = h5read(data.h5file,'/Image',[1 indx 1],[data.sz(1) 1 data.sz(3)]);
    set(data.hImages(3), 'CData', squeeze(Iyz));
    if (isfield(data,'Vec') && ...
            (pos >= data.Vec.x(1)) && (pos <= data.Vec.x(end))),
        vindx = (indx - 1)/(data.sz(2)-1) * (size(data.Vec.u,2)-1) + 1;
        
        delete(data.Vec.h{3});
        
        hold on;
        data.Vec.h{3} = quiverc(data.Vec.z+data.Vec.yzOffset,data.Vec.y, ...
            data.Vec.w(:,vindx,:),data.Vec.v(:,vindx,:), ...
            'AbsScale',data.Vec.scale,'g');
        hold off;
    end;
end

if ~isempty(indz)
    data.indz = indz;
    pos1 = data.Z{1}(1) + (indz - 1)/(data.sz(3)-1) * diff(data.Z{1});
    pos2 = data.Z{2}(1) + (indz - 1)/(data.sz(3)-1) * diff(data.Z{2});
    
    % adjust the line position
    set(data.hLines(3,1),'YData',[pos1 pos1]);
    set(data.hLines(3,2),'XData',[pos2 pos2]);
    
    % and get the new plane
    Ixy = h5read(data.h5file,'/Image',[1 1 indz],[data.sz(1) data.sz(2) 1]);
    set(data.hImages(1), 'CData', Ixy);
    
    if (isfield(data,'Vec') && ...
            (pos >= data.Vec.z(1)) && (pos < data.Vec.z(end))),
        vindz = (indz - 1)/(data.sz(3)-1) * (size(data.Vec.u,3)-1) + 1;
        
        delete(data.Vec.h{1});
        
        hold on;
        data.Vec.h{1} = quiverc(data.Vec.x,data.Vec.y, ...
            data.Vec.u(:,:,vindz),data.Vec.v(:,:,vindz), ...
            'AbsScale',data.Vec.scale,'y');
        hold off;
    end;
end;

% -------------------------------------------------------------
function orthoLineDown(obj,evendata,ln)
% We clicked on a line, so now set up the button motion handler to track
% where we drag it to

data = guidata(obj);

data.dragAxes = get(obj,'Parent');
set(data.Figure,'WindowButtonMotionFcn',{@orthoLineMove,ln},...
                'WindowButtonUpFcn',{@orthoLineUp,ln});

guidata(obj,data);

% -------------------------------------------------------------
function orthoLineUp(obj,evendata,ln)
% The important mouse function.  Figure out where we clicked and adjust the
% planes appropriately.  The parameter ln indicates which line we clicked
% on.
data = guidata(obj);

c = get(data.dragAxes, 'CurrentPoint');

set(data.Figure,'WindowButtonMotionFcn','', ...
    'WindowButtonUpFcn','');

indx = [];
indy = [];
indz = [];
switch ln,
    case 1,                                % xz plane
        pos = c(1,2);
        if ((pos >= data.Y(1)) && (pos <= data.Y(2))),
            indy = round((pos - data.Y(1))/diff(data.Y)*(data.sz(1)-1))+1;
        end;
        
    case 2,
        pos = c(1,1);                         % yz plane
        if ((pos >= data.X(1)) && (pos <= data.X(2))),
            indx = round((pos - data.X(1))/diff(data.X)*(data.sz(2)-1))+1;
        end;
        
    case 3,                                % xy plane
        % check which of the two lines we clicked (in the xz or yz views)
        if ((c(1,2) >= data.Z{1}(1)) && (c(1,2) <= data.Z{1}(2))),
            pos = c(1,2);
            pos = pos-data.Z{1}(1);
        elseif ((c(1,1) >= data.Z{2}(1)) && (c(1,1) <= data.Z{2}(2))),
            pos = c(1,1);
            pos = pos-data.Z{2}(1);
        else
            return;
        end;
        
        % and get the new plane
        indz = round(pos/diff(data.Z{1})*(data.sz(3)-1))+1;
end;

data = orthoUpdate(data,indx,indy,indz);

guidata(obj,data);

% -------------------------------------------------------------
function orthoLineMove(obj,evendata,ln)
% once the mouse comes up, remove the button handler functions

data = guidata(obj);
c = get(data.dragAxes, 'CurrentPoint');

switch ln,
    case 1,                                % xz plane
        pos = c(1,2);
        if ((pos >= data.Y(1)) && (pos <= data.Y(2))),
            set(data.hLines(ln,1),'YData',[pos pos]);
        end
        
    case 2,                             % yz plane
        pos = c(1,1);                         % yz plane
        if ((pos >= data.X(1)) && (pos <= data.X(2))),
            set(data.hLines(ln,1),'XData',[pos pos]);
        end
        
    case 3,                                 %xy plane
        % check which of the two lines we clicked (in the xz or yz views)
        if ((c(1,2) >= data.Z{1}(1)) && (c(1,2) <= data.Z{1}(2))),
            pos = c(1,2);
            pos = pos-data.Z{1}(1);
        elseif ((c(1,1) >= data.Z{2}(1)) && (c(1,1) <= data.Z{2}(2))),
            pos = c(1,1);
            pos = pos-data.Z{2}(1);
        else
            return;
        end;
        % adjust the line position
        set(data.hLines(ln,1),'YData',[pos pos]+data.Z{1}(1));
        set(data.hLines(ln,2),'XData',[pos pos]+data.Z{2}(1));
end

guidata(obj,data);

% -------------------------------------------------------------
function orthoKeyPress(obj,eventdata)

data = guidata(obj);

indx = [];
indy = [];
indz = [];
switch eventdata.Key
    case 'u'
        if data.indz < data.sz(3)
            indz = data.indz + 1;
        end
    case 'd'
        if data.indz > 1
            indz = data.indz - 1;
        end
        
    case 'f'
        if data.indx < data.sz(2)
            indx = data.indx + 1;
        end
    case 'b'
        if data.indx > 1
            indx = data.indx - 1;
        end
        
    case 'l'
        if data.indy < data.sz(1)
            indy = data.indy + 1;
        end
    case 'r'
        if data.indy > 1
            indy = data.indy - 1;
        end
end

if ~isempty(indx) || ~isempty(indy) || ~isempty(indz)
    data = orthoUpdate(data,indx,indy,indz);
end

guidata(obj,data);

        
        


