function h = orthoplot(a,b,c, I, varargin)
% function h = orthoplot(a,b,r, I, vx,vy,vz, u,v,w)
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

opt.order = 'ijk';

if ((nargin >= 10) && ...
    isnumeric(varargin{1}) && isnumeric(varargin{2}) && isnumeric(varargin{3}) && ...
    isnumeric(varargin{4}) && isnumeric(varargin{5}) && isnumeric(varargin{6}))
    [vx,vy,vz, u,v,w] = deal(varargin{1:6});
    varargin = varargin(11:end);
    p = 11;
else
    p = 5;
end

opt = parsevarargin(opt,varargin, p);

if (nargin == 1),
    I = a;
    
    switch opt.order
        case 'xyz'
            a = 1:size(I,2);
            b = 1:size(I,1);
            c = 1:size(I,3);
        case 'ijk'
            a = 1:size(I,1);
            b = 1:size(I,2);
            c = 1:size(I,3);
            I = permute(I,[2 1 3]);
    end
elseif (ndims(a) == 3),
    switch opt.order
        case 'xyz'
            % do nothing
        case 'ijk'
            if all(size(a) ~= 1)
                a = permute(a,[2 1 3]);
                b = permute(b,[2 1 3]);
                c = permute(c,[2 1 3]);
                I = permute(I,[2 1 3]);
            end
    end
    if (any(flatten(diff(a,[],1)) ~= 0) || ...
        any(flatten(diff(a,[],3)) ~= 0) || ...
        any(flatten(diff(b,[],2)) ~= 0) || ...
        any(flatten(diff(b,[],3)) ~= 0) || ...
        any(flatten(diff(c,[],1)) ~= 0) || ...
        any(flatten(diff(c,[],2)) ~= 0)),
        error('x,y and z must be vectors or plaid.');
    end;
    % and turn the matrices into vectors
    a = a(1,:,1);
    b = b(:,1,1);
    c = flatten(c(1,1,:));
end;

if (nargin >= 10) && strcmp(opt.order,'ijk')
    vx1 = permute(vy,[2 1 3]);
    vy1 = permute(vx,[2 1 3]);
    vx = vx1;
    vy = vy1;
    vz = permute(vz,[2 1 3]);

    u1 = permute(v,[2 1 3]);
    v1 = permute(u,[2 1 3]);
    u = u1;
    v = v1;
    w = permute(w,[2 1 3]);
end

if (ndims(I) ~= 3),
  error('Image data must be three dimensional.');
end;

% define the x and y positions of each orthogonal view
% since we're all in one axis, the xz and yz views have to be offset
% down and right, respectively.
xyx = a([1 end]);
xyy = b([1 end]);
xzx = a([1 end]);
xzy = [0 range(c)] + b(end);
yzx = [0 range(c)] + a(end);
yzy = b([1 end]);

clim = [min(I(:)) max(I(:))];

% get the middle index from each coordinate array
% TODO: in principle, we shouldn't assume even spacing, as we do here
mid = ceil([length(a)/2 length(b)/2 length(c)/2]);

% plot the images
himxy = imagesc(xyx,xyy, I(:,:,mid(3)));
hold on;
himxz = imagesc(xzx,xzy, squeeze(I(mid(2),:,:))');
himyz = imagesc(yzx,yzy, squeeze(I(:,mid(1),:)));

% yes, the "tight" needs to be there twice.  Don't ask me why, but it does
% what I want it to that way
axis tight equal tight off;
caxis(clim);
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
    [q,mx] = min(abs(vx-(a(end)+a(1))/2));
    [q,my] = min(abs(vy-(b(end)+b(1))/2));
    [q,mz] = min(abs(vz-(c(end)+c(1))/2));

    mid = {mx my mz};

    if all(vx == vx(1)) || all(vy == vy(1)) || all(vz == vz(1))
        mid = {1:length(vx), 1:length(vy), 1:length(vz)};
    end
    
    % offset for the xz and yz axes
    xzOffset = max(b) - min(c);
    yzOffset = max(a) - min(c);

    % calculate a 3D scaling value
    mag = sqrt(u.^2 + v.^2 + w.^2);
    if (numel(vx) > 1),
        d = nanmean(cat(1,diff(vx(:)),diff(vy(:)),diff(vz(:))));
    else
        d = min(size(I))/2;
    end;
    if (d == 0)
        d = min([range(a) range(b) range(c)]);
    end
    vscale = d/prctile(mag(:),95);

    % and draw the vectors
    hvecxy = quiverc(vx,vy, ...
                     u(:,:,mid{3}),v(:,:,mid{3}),'y',...
                     'AbsScale',vscale);
    hvecxz = quiverc(vx,vz+xzOffset, ...
                     squeeze(u(mid{2},:,:)),...
                     squeeze(w(mid{2},:,:)),'r','AbsScale',vscale);
    hvecyz = quiverc(vz+yzOffset,vy,...
                     squeeze(w(:,mid{1},:)),...
                     squeeze(v(:,mid{1},:)),'g','AbsScale',vscale);

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
set(data.Figure, 'DoubleBuffer', 'on');

data.X = xyx;
data.Y = xyy;
data.Z = {xzy yzx};
data.I = I;

data.hImages = [himxy himxz himyz];
data.hLines = hline;
guidata(data.Figure, data);

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
function orthoLineMove(obj,evendata,ln)
% The important mouse function.  Figure out where we clicked and adjust the
% planes appropriately.  The parameter ln indicates which line we clicked
% on.
data = guidata(obj);
c = get(data.dragAxes, 'CurrentPoint');

switch ln,
 case 1,                                % xz plane
  pos = c(1,2);
  if ((pos >= data.Y(1)) && (pos <= data.Y(2))),
      set(data.hLines(ln,1),'YData',[pos pos]);
      
      ind = round((pos - data.Y(1))/diff(data.Y)*(size(data.I,1)-1))+1;
      set(data.hImages(2), 'CData', squeeze(data.I(ind,:,:))');
  end;
  if (isfield(data,'Vec') && ...
      (pos >= data.Vec.y(1)) && (pos <= data.Vec.y(end))),
      ind = round((pos - data.Vec.y(1))/range(data.Vec.y) * ...
                  (size(data.Vec.u,1)-1))+1;
      
      delete(data.Vec.h{2});
      
      hold on;
      data.Vec.h{2} = quiverc(data.Vec.x,data.Vec.z+data.Vec.xzOffset, ...
                              data.Vec.u(ind,:,:),data.Vec.w(ind,:,:), ...
                              'AbsScale',data.Vec.scale,'r');
      hold off;
  end;

 case 2,
  pos = c(1,1);                         % yz plane
  if ((pos >= data.X(1)) && (pos <= data.X(2))),
      set(data.hLines(ln,1),'XData',[pos pos]);
    
      ind = round((pos - data.X(1))/diff(data.X)*(size(data.I,2)-1))+1;
      set(data.hImages(3), 'CData', squeeze(data.I(:,ind,:)));
  end;
  if (isfield(data,'Vec') && ...
      (pos >= data.Vec.x(1)) && (pos <= data.Vec.x(end))),
      ind = round((pos - data.Vec.x(1))/range(data.Vec.x) * ...
                  (size(data.Vec.u,2)-1))+1;
      
      delete(data.Vec.h{3});
      
      hold on;
      data.Vec.h{3} = quiverc(data.Vec.z+data.Vec.yzOffset,data.Vec.y, ...
                              data.Vec.w(:,ind,:),data.Vec.v(:,ind,:), ...
                              'AbsScale',data.Vec.scale,'g');
      hold off;
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

  % adjust the line position
  set(data.hLines(ln,1),'YData',[pos pos]+data.Z{1}(1));
  set(data.hLines(ln,2),'XData',[pos pos]+data.Z{2}(1));

  % and get the new plane
  ind = round(pos/diff(data.Z{1})*(size(data.I,3)-1))+1;
  set(data.hImages(1), 'CData', squeeze(data.I(:,:,ind)));

  if (isfield(data,'Vec') && ...
      (pos >= data.Vec.z(1)) && (pos < data.Vec.z(end))),
      ind = round((pos - data.Vec.z(1))/range(data.Vec.z) * ...
                  (size(data.Vec.u,3)-1))+1;
      
      delete(data.Vec.h{1});
      
      hold on;
      data.Vec.h{1} = quiverc(data.Vec.x,data.Vec.y, ...
                              data.Vec.u(:,:,ind),data.Vec.v(:,:,ind), ...
                              'AbsScale',data.Vec.scale,'y');
      hold off;
  end;
end;

guidata(obj,data);

% -------------------------------------------------------------
function orthoLineUp(obj,evendata,ln)
% once the mouse comes up, remove the button handler functions

data = guidata(obj);

set(data.Figure,'WindowButtonMotionFcn','', ...
                'WindowButtonUpFcn','');

guidata(obj,data);

