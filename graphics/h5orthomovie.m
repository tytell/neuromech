function h5orthomovie(h5file,outfile, varargin)

opt.fps = 5;
opt.planes = 'xyz';
opt.outsize = 640;
opt.range = [1 Inf; 1 Inf; 1 Inf];
opt.datasetname = '/Image';

opt = parsevarargin(opt,varargin, 2);

info = h5info(h5file,opt.datasetname);
sz = info.Dataspace.Size;
fprintf('Total size: %d y, %d x, %d z\n', sz);

scale = 1;
if ~isempty(info.Attributes)
    attrnames = {info.Attributes.Name};
    if ismember('Scale',attrnames)
        scale = h5readatt(h5file,opt.datasetname,'Scale');
    end
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
mid = ceil([length(y)/2 length(x)/2 length(z)/2]);

if any(opt.planes == 'x')
    indx = opt.range(2,1);
    if (opt.range(2,2) > sz(2))
        opt.range(2,2) = sz(2);
    end
    dx = 1;
else
    indx = mid(2);
    dx = 0;
end
if any(opt.planes == 'y')
    indy = opt.range(1,1);
    if (opt.range(1,2) > sz(1))
        opt.range(1,2) = sz(1);
    end
    dy = 1;
else
    indy = mid(1);
    dy = 0;
end
if any(opt.planes == 'z')
    indz = opt.range(3,1);
    if (opt.range(3,2) > sz(3))
        opt.range(3,2) = sz(3);
    end
    dz = 1;
else
    indz = mid(3);
    dz = 0;
end

Ixy = h5read(h5file,opt.datasetname,[1 1 indz],[sz(1) sz(2) 1]);
Ixz = h5read(h5file,opt.datasetname,[indy 1 1],[1 sz(2) sz(3)]);
Iyz = h5read(h5file,opt.datasetname,[1 indx 1],[sz(1) 1 sz(3)]);

% plot the images
hax = axes('Position',[0.02 0.02 0.96 0.96]);

himxy = imagesc(xyx,xyy, Ixy, 'Parent',hax);
hold on;
himxz = imagesc(xzx,xzy, squeeze(Ixz)', 'Parent',hax);
himyz = imagesc(yzx,yzy, squeeze(Iyz), 'Parent',hax);

% yes, the "tight" needs to be there twice.  Don't ask me why, but it does
% what I want it to that way
axis(hax,'tight','equal','tight','off');
colormap gray;

% box in the different views
plot(xyx([1 2 2 1 1]),xyy([1 1 2 2 1]),'y-',...
     xzx([1 2 2 1 1]),xzy([1 1 2 2 1]),'r-',...
     yzx([1 2 2 1 1]),yzy([1 1 2 2 1]),'g-',...
     'LineWidth',1.5);

% and draw the lines that indicate where the planes are
hline = repmat(-1,[3 2]);
hline(1,1) = plot([xyx(1) yzx(2)],mean(xyy)*[1 1],'r');
hline(2,1) = plot(mean(xyx)*[1 1],[xyy(1) xzy(2)],'g');
hline(3,1) = plot(xyx,mean(xzy)*[1 1],'y');
hline(3,2) = plot(mean(yzx)*[1 1],yzy,'y');

hold off;

% set up the guidata for the button handler routines
data.Figure = gcf;
set(data.Figure,'Color','w', 'Units','pixels','WindowStyle','normal');

pos = get(data.Figure,'Position');
width = pos(3);
height = pos(4);
fac = opt.outsize / width;
width2 = width*fac;
height2 = height*fac;

pos2 = [pos(1) pos(2)+height-height2 width2 height2];
set(data.Figure, 'Position',pos2);

data.X = xyx;
data.Y = xyy;
data.Z = {xzy yzx};
data.indx = indx;
data.indy = indy;
data.indz = indz;
data.h5file = h5file;
data.sz = sz;

data.hImages = [himxy himxz himyz];
data.hLines = hline;

if ~isempty(outfile)
    vid = VideoWriter(outfile);
    set(vid, 'FrameRate',opt.fps);
    open(vid);
end

N = max(diff(opt.range,[],2));
timedWaitBar(0,'Plotting images...');
for i = 1:N
    data.indx = data.indx + dx;
    if (data.indx == opt.range(2,1))
        dx = 1;
    elseif (data.indx == opt.range(2,2))
        dx = -1;
    end
    data.indy = data.indy + dy;
    if (data.indy == opt.range(1,1))
        dy = 1;
    elseif (data.indy == opt.range(1,2))
        dy = -1;
    end
    data.indz = data.indz + dz;
    if (data.indz == opt.range(3,1))
        dz = 1;
    elseif (data.indz == opt.range(3,2))
        dz = -1;
    end
    
    data = orthoUpdate(data, data.indx,data.indy,data.indz);
    drawnow;
    
    if ~isempty(outfile)
        frame = getframe;
        writeVideo(vid,frame);
    end
    timedWaitBar(i/N);
end
if ~isempty(outfile)
    close(vid);
end
timedWaitBar(1);

% -------------------------------------------------------------
function data = orthoUpdate(data, indx,indy,indz)

if ~isempty(indy)
    data.indy = indy;
    
    pos = data.Y(1) + (indy - 1)/(data.sz(1)-1) * diff(data.Y);
    set(data.hLines(1,1),'YData',[pos pos]);
    
    Ixz = h5read(data.h5file,opt.datasetname,[indy 1 1],[1 data.sz(2) data.sz(3)]);
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
    
    Iyz = h5read(data.h5file,opt.datasetname,[1 indx 1],[data.sz(1) 1 data.sz(3)]);
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
    Ixy = h5read(data.h5file,opt.datasetname,[1 1 indz],[data.sz(1) data.sz(2) 1]);
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


