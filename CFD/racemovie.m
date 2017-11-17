function racemovie(files,varargin)
% function racemovie(files,varargin)
% Builds up a movie of one or more swimmers swimming together.  Used to
% generate the race movie that I tend to show in presentations.
%
% Options -
%   'showall' - Show all frames for all movies, even if some are shorter.
%     Otherwise stops when the shortest movie ends.
%   'outputfile' - Generates an AVI in this file.  If empty, do not
%     generate the AVI.


opt.showall = false;
opt.outputfile = '';
opt.fps = 5;
opt.delay = 0;
opt.videosize = [1024 512];
opt.colors = [];
opt = parsevarargin(opt,varargin,1);

t = [];
ox = [];
oy = [];
for i = 1:length(files),
    F = load(files{i},'t','xm','ym','xl','yl','xr','yr');
    
    if (getvar('-file',files{i},'good')),
        F.t = F.t(good);
        F.xm = F.xm(:,good);
        F.ym = F.ym(:,good);
        F.xl = F.xl(:,good);
        F.yl = F.yl(:,good);
        F.xr = F.xr(:,good);
        F.yr = F.yr(:,good);
    end;
    
    [pn,fn] = fileparts(files{i});
    load(fullfile(pn,[fn '_analysis.mat']),'swimvecx','swimvecy');

    npt = size(F.xm,1);
    nfr = size(F.xm,2);
    
    vx = -1; %mean(swimvecx(t > 2));
    vy = 0; %mean(swimvecy(t > 2));
    
    ox1 = cat(1,F.xl,F.xr(end:-1:1,:));
    oy1 = cat(1,F.yl,F.yr(end:-1:1,:));
    ox2 = -ox1.*vx - oy1.*vy;
    oy2 = -ox1.*vy + oy1.*vx + 4*i;
    
    t = catuneven(3,t,F.t);
    ox = catuneven(3,ox,ox2);
    oy = catuneven(3,oy,oy2);
end;

dt = t(:,2,:) - t(:,1,:);
[~,ind] = max(dt);
t0 = t(:,:,ind);
good0 = isfinite(t0);
for i = 1:length(files),
    good = isfinite(t(:,:,i));
    for j = 1:size(ox,1),
        ox(j,good0,i) = interp1(t(:,good,i),ox(j,good,i), t0(:,good0));
        oy(j,good0,i) = interp1(t(:,good,i),oy(j,good,i), t0(:,good0));
    end;
    ox(:,~good0,i) = NaN;
    oy(:,~good0,i) = NaN;
end;
ox = ox(:,good0,:);
oy = oy(:,good0,:);

if (~opt.showall),
    good = all(isfinite(ox(1,:,:)),3);
    ox = ox(:,good,:);
    oy = oy(:,good,:);
end;

axlim = [min(ox(:)) max(ox(:)) min(oy(:)) max(oy(:))];
nfr = size(ox,2);

ox = permute(ox,[1 3 2]);
oy = permute(oy,[1 3 2]);

if (~isempty(opt.outputfile)),
    aviobj = VideoWriter(opt.outputfile);
    aviobj.FrameRate = opt.fps;
    open(aviobj);
end;

nr = size(ox,2);

fig = figure('WindowStyle','normal', 'Units','pixels');
pos = get(fig, 'Position');

t = pos(2) + pos(4);
b = t - opt.videosize(2);
pos = [pos(1) b opt.videosize];

set(fig, 'Position', pos, 'Color','w');

clf;
ax = axes('Position',[0 0 1 1]);
hold(ax,'on');

h = zeros(1,size(ox,2));
for j = 1:nr,
    h(j) = fill(ox(:,j,1),oy(:,j,1),j, 'EdgeColor','none', 'Parent',ax);
    if ~isempty(opt.colors)
        set(h(j), 'FaceColor', opt.colors(j,:));
    end
end;
hold(ax,'off');

for i = 1:nfr,
    for j = 1:nr,
        set(h(j),'XData',ox(:,j,i), 'YData',oy(:,j,i));
    end;
    
    axis(ax,'equal','off');
    axis(ax, axlim);
    
    drawnow;
    
    if (~isempty(opt.outputfile)),
        F = getframe(fig);
        sz = size(F.cdata);
        sz = floor(sz(1:2)/4)*4;
        F.cdata = F.cdata(1:sz(1),1:sz(2),:);
        
        writeVideo(aviobj,F);
    end;
        
    if (opt.delay > 0),
        pause(opt.delay);
    end;
end;

if (~isempty(opt.outputfile)),
    close(aviobj);
end;
