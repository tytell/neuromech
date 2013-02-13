function [mx,my, lx,ly, rx,ry] = getMidlineAvi(aviname, frames, hx,hy, ...
                                              tx,ty, n, width, bright, ...
                                              quiet, maxang)

if (nargin < 11),
    maxang = [];
    if (nargin < 10),
        quiet = 0;
        if (nargin < 9),
            bright = 1;
        end;
    end;
end;

fig = gcf;
set(fig, 'DoubleBuffer', 'on');

set(fig, 'CurrentCharacter', '0');

k = find(isfinite(hx) & isfinite(hy) & isfinite(tx) & isfinite(ty));

mx = repmat(NaN,[n length(hx)]);
my = repmat(NaN,[n length(hx)]);
lx = repmat(NaN,[n length(hx)]);
ly = repmat(NaN,[n length(hx)]);
rx = repmat(NaN,[n length(hx)]);
ry = repmat(NaN,[n length(hx)]);

if (quiet < 2),
    timedWaitBar(0, 'Finding midlines...');
end;

for ii = 1:length(k),
    i = k(ii);
    
    I = im2double(frame2im(aviread(aviname, frames(i))));
    if (~bright),
        I = 1-I;
    end;
    
    [mx1,my1, lx1,ly1, rx1,ry1] = getMidline(I, hx(i),hy(i), ...
                                                tx(i),ty(i), n,width, maxang);
    
    if (quiet == 0),
        imshow(I,'n');
        addplot(mx1,my1,'b.-', lx1,ly1,'g.-', rx1,ry1,'r.-');
        drawnow;
        
        if (lower(get(fig,'CurrentCharacter')) == 'q'),
            break;
        end;
    end;
    
    mx(:,i) = mx1;
    my(:,i) = my1;
    lx(:,i) = lx1;
    ly(:,i) = ly1;
    rx(:,i) = rx1;
    ry(:,i) = ry1;
    
    if (quiet < 2),
        if (~timedWaitBar(ii/length(k))),
            break;
        end;
    end;
end;
