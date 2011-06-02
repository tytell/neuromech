function [mx,my,good] = getMidline2(avi,frames, mx1,my1, hx,hy, npt, ...
                                        width,len,maxang,varargin)

global GMdebug;

mincorrel = 0.4;
peakheight = 0.95;
maxpeaksz = 0.4;

info = aviinfo(avi);
if (isempty(frames)),
    frames = 1:info.NumFrames;
elseif (numel(frames) == 1),
    frames = 1:frames:info.NumFrames;
end;
nfr = info.NumFrames;

opts = varargin;

invert = logical(0);
quiet = logical(0);
isvert = logical(1);
GMdebug = 0;
facedir = 1;
reduce = 1;
i = 1;
while (i <= length(opts)),
    switch lower(opts{i}),
     case 'invert',
      invert = 1;
     case 'quiet',
      quiet = 1;
     case {'horizontal','horiz'},
      isvert = 0;
     case {'facedown','faceright'},
      facedir = 1;
     case {'faceup','faceleft'},
      facedir = -1;
     case 'debug',
      if ((i < length(opts)) & isnumeric(opts{i+1})),
          GMdebug = opts{i+1};
          i = i+1;
      else
          GMdebug = 1;
      end;
     case 'reduce',
      reduce = opts{i+1};
      i = i+1;
     otherwise,
      error('Unrecognized option %s',opts{i});
    end;
    i = i+1;
end;

if (reduce < 1),
    mx1 = mx1.*reduce;
    my1 = my1.*reduce;
    hx = hx.*reduce;
    hy = hy.*reduce;
    len = len.*reduce;
end;

width = width*len;
width = spline(linspace(0,len,length(width)),width, ...
               linspace(0,len,npt));

seglen = len/(npt-1);
s0 = (0:seglen:len)';
perpd = 2*tan(maxang);

if (isempty(mx1)),
    I = frame2im(aviread(avi,frames(1)));
    if (reduce < 1),
        I = imresize(I,reduce,'bilinear');
    end;

    imshow(I,'n');
    set(gcf,'Name','Click any number of points along first midline');
    
    [mx1,my1] = ginputb;
    mx1 = mx1';
    my1 = my1';
end;

if ((length(mx1) ~= npt) | ...
    any((sqrt(diff(mx1).^2 + diff(my1).^2)-seglen)/seglen > 0.05)),

    mxclick = mx1;
    myclick = my1;

    sclick = [0; cumsum(sqrt(diff(mxclick).^2 + diff(myclick).^2))];

    if (abs(sclick(end)-len)/len > 0.05),
        warning(['Clicked length is more than 5% different that input ' ...
                 'length']);
    end;

    sp = spaps(sclick,[mxclick myclick]',0.25);
    k = find(s0 <= sclick(end));
    mxy = fnval(sp,s0(k));
    mx1(k) = mxy(1,:)';
    my1(k) = mxy(2,:)';

    ang1 = repmat(NaN,[npt-1 1]);
    ang1(k(2:end)) = atan2(diff(my1(k)),diff(mx1(k)));
    ang1(1) = ang1(2);

    if (length(k) < npt),
        newpt = k(end)+1:npt;
        oldpt = k(end-4:end);

        ff = repmat(frames(1),[length(oldpt) 1]);
        ang1(newpt) = ...
            extrap2(ff,s0(oldpt),ang1(oldpt),2, ...
                          repmat(frames(1),[length(newpt),1]),s0(newpt));
    end;
end;

ang = repmat(NaN,[npt nfr]); 
ang(:,frames(1)) = ang1;

good = zeros(npt,nfr);

ff = repmat(NaN,[npt nfr]);
ff(:,frames) = repmat(frames,[npt 1]);

for j = 1:length(frames),
    fr = frames(j);

    I = im2double(frame2im(aviread(avi,fr)));
    if (size(I,3) > 1),
        I = I(:,:,1);
    end;
    if (invert),
        I = 1-I;
    end;
    if (reduce < 1),
        I = imresize(I,reduce,'bilinear');
    end;

    mx1(1) = hx(fr);
    mx1(2:npt) = NaN;
    my1(1) = hy(fr);
    my1(2:npt) = NaN;
    good(1,fr) = 1;
    lastgood = 1;
    for i = 2:npt,
        dxds = cos(ang1(i-1));
        dyds = sin(ang1(i-1));

        mx1(i) = mx1(i-1) + seglen*dxds;
        my1(i) = my1(i-1) + seglen*dyds;

        w1 = abs(1.5*dyds*width(i)*(i-lastgood)) + abs(0.5*seglen*dxds);
        h1 = abs(1.5*dxds*width(i)*(i-lastgood)) + abs(0.5*seglen*dyds);

        rgnx = round(mx1(i) + [-w1 w1]/2);
        rgny = round(my1(i) + [-h1 h1]/2);
        if ((rgnx(1) < 1) | (rgnx(2) > size(I,2)) | ...
            (rgny(1) < 1) | (rgny(2) > size(I,1))),
            warning('A portion of the fish has gone off the image.');
            break;
        end;

        Irgn = I(rgny(1):rgny(2),rgnx(1):rgnx(2));

        angguess = nanmean(ang1([i-1 i]));
        [dpos,newang,maxcorr,peaksz] = ...
            correl(Irgn,width(i), angguess,maxang);
        
        if (GMdebug > 1),
            clf;
            imshow(I,'n');
            addplot(mx1,my1,'r.-',...
                    mx1(i) + seglen*cos(ang1(i)),...
                    my1(i) + seglen*sin(ang1(i)),'g.')
            drawnow;
            
            if (i == 17),
                GMdebug = 4;
            end;
        end;
        if ((maxcorr >= mincorrel) & (peaksz < maxpeaksz)),
            dxds2 = cos(newang);
            dyds2 = sin(newang);

            newx = mx1(i) - dpos*dyds2;
            newy = my1(i) + dpos*dxds2;
        
            ang1(i-1) = atan2(newy - my1(i-1), newx - mx1(i-1));
            ang1(i) = newang;
        
            mx1(i) = mx1(i-1) + seglen*cos(ang1(i-1));
            my1(i) = my1(i-1) + seglen*sin(ang1(i-1));

            good(i,fr) = 1;
            lastgood = i;
        else
            good(i,fr) = 0;
            ang(:,fr) = ang1;

            if (j > 3),
                indfr = j-3:j+2;
            else
                indfr = j;
            end;
            
            goodpt = find(good(:,fr));
            if (length(goodpt) > 4),
                indpos = goodpt(end-3:end);
            else
                indpos = goodpt;
            end;
            
            goodang = ang(indpos,frames(indfr));
            goodang(~good(indpos,frames(indfr))) = NaN;
            newang = extrap2(frames(indfr),s0(indpos),...
                           goodang,2, fr,s0(i));
            ang1(i) = newang;
        end;
    end;

    k = find(good(:,fr));
    bad = find(~good(:,fr));
    if (~isempty(bad) & any(bad < k(end))),
        bad = bad(bad < k(end));
        ang(bad,fr) = spline(s0(k),ang(k,fr),s0(bad));
    end;

    if (~quiet),
        imshow(I,'n');
        drawx = mx1;
        drawy = my1;
        drawx(~good(:,fr)) = NaN;
        drawy(~good(:,fr)) = NaN;
        addplot(drawx,drawy,'r.-',mx1(~good(:,fr)),my1(~good(:,fr)),'g.');
        drawnow;

        set(gcf,'Name',sprintf('Frame %d/%d (%d%%)',fr,nfr, ...
                               round(j/length(frames)*100)));
    end;
end;

mx = repmat(NaN,[npt length(frames)]);
my = repmat(NaN,[npt length(frames)]);



% ----------------------------------------------------------------
function [dpos,newang,maxcorr,peaksize] = correl(I,segwidth, ang0,maxdang)

global GMdebug;

d1 = min(size(I));
d2 = max(size(I));
% angular distance corresponding to one pixel at the edge of the image
angstep = 2/d1;                          

w2 = segwidth/2;

% only take at most 40 steps (from -maxdang to maxdang)
if (maxdang/angstep > 20),
    angstep = maxdang/20;
end;

ang = -maxdang:angstep:maxdang;
% center around 0
ang = ang0 + ang + (maxdang-ang(end))/2;
nang = length(ang);

[x,y] = meshgrid((1:size(I,2)) - (size(I,2)/2),...
                (1:size(I,1)) - (size(I,1)/2));

C = cos(ang);
S = sin(ang);

maxdisp = floor(max(cat(2,abs(-S*x(1,1)+C*y(1,1)),...
                        abs(-S*x(1,end)+C*y(1,end)),...
                        abs(-S*x(end,1)+C*y(end,1)),...
                        abs(-S*x(end,end)+C*y(end,end)))));
dy = (-maxdisp+floor(w2)):(maxdisp-floor(w2));
ndisp = length(dy);

if (prod([size(I) nang ndisp]) > 1e7),
    error('Correlation region is too large.  Please reduce the image.');
end;

template = zeros(size(I,1),size(I,2),nang,ndisp);

Call = repmat(reshape(C,[1 1 nang]),[size(x) 1]);
Sall = repmat(reshape(S,[1 1 nang]),[size(x) 1]);
yrot = repmat(-Sall.*repmat(x,[1 1 nang]) + Call.*repmat(y,[1 1 nang]), ...
              [1 1 1 ndisp]);

dyall = repmat(reshape(dy,[1 1 1 ndisp]),[size(x) nang 1]);

% ymax = floor(max(abs(flatten(yrot,[1 2]))));
% ymax = shiftdim(ymax,-1);
% ymax = repmat(ymax,[size(x) 1 1]) - floor(w2);

% dy(abs(dy) > ymax) = NaN;

template(abs(yrot-dyall) <= w2) = 1;
N = sum(flatten(template,[1 2]));
N(N == 0) = NaN;

Iall = repmat(I,[1 1 nang ndisp]);

corr = sum(flatten(Iall.*template, [1 2]))./N;
corr = squeeze(corr);


% for i = 1:length(ang),
%     yrot = -S(i)*x + C(i)*y;

%     ymax = floor(max(abs(yrot(:))));
%     dy1 = (-ymax+floor(w2)):(ymax-floor(w2));
%     corr1 = zeros(size(dy1));
%     maxcorr1 = 0;

%     a = find(dy == -ymax+floor(w2)) - 1;

%     for j = 1:length(dy1),
%         Ic = zeros(size(I));
%         Ic(abs(yrot-dy1(j)) <= w2) = 1;

%         corr1(j) = sum(sum(Ic.*I))/sum(Ic(:));
%     end;

%     corr(i,(1:length(dy1))+a) = corr1;
% end;

[maxcorr,k] = max(corr(:));
[i,j] = ind2sub(size(corr),k);
peaksize = sum(corr(:) > 0.95*maxcorr)/numel(corr);

dpos = dy(j);
newang = ang(i);

if (GMdebug > 3),
    subplot(2,1,1);
    imagesc(dy,ang*180/pi,corr);
    colormap default;
    addplot(dpos,newang*180/pi,'go');

    subplot(2,1,2);
    yrot = -S(i)*x + C(i)*y;
    Ic = zeros(size(I));
    Ic(abs(yrot-dy(j)) <= w2) = 1;
    Idbg(:,:,1) = I;
    Idbg(:,:,2) = I.*(1-Ic);
    Idbg(:,:,3) = I.*(1-Ic);
    imshow(Idbg,'n');
    title('True max');
    axis on;
    pause;
elseif (GMdebug > 2),
    imshow(x(1,:),y(:,1),I,'n');

    ctrx = -dpos*sin(newang);
    ctry = dpos*cos(newang);
    addplot([0 cos(ang0)]*d1/3,[0 sin(ang0)]*d1/3,'o-',...
            [0 cos(newang)]*d1/3 + ctrx, [0 sin(newang)]*d1/3 + ctry,'o-');
    drawnow;
end;
