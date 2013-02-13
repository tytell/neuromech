function streakImage(in,out,n,frames,crop)

a = 3;

info = aviinfo(in);
if (nargin < 5),
    crop = [1 info.Width 1 info.Height];
    if (nargin < 4),
        frames = [];
    end;
end;
if (isempty(frames)),
    frames = 1:info.NumFrames-n+1;
end;
if (frames(end) > info.NumFrames-n+1),
    warning(['Frame list is longer than AVI.  Reducing number of ' ...
             'frames.']);
    k = find(frames <= info.NumFrames-n+1);
    frames = frames(k);
end;

rcrop = crop(3):crop(4);
ccrop = crop(1):crop(2);

w = length(ccrop);
h = length(rcrop);

fr0 = aviread(in,(1:n-1)+frames(1));
for i = 1:n-1,
    fr0(i).cdata = fr0(i).cdata(rcrop,ccrop);
end;
fr(1) = fr0(1);
fr(2:n) = fr0;

emov = avifile(out,'Compression','none');

fr2.cdata = uint8(zeros(h,w));
fr2.colormap = fr(1).colormap;

timedWaitBar(0,'Creating streaks');
for ii = 1:length(frames),
    i = frames(ii);

    fr(1:n-1) = fr(2:n);
    fr(n) = aviread(in,i+n-1);
    fr(n).cdata = fr(n).cdata(rcrop,ccrop);

    I = a*sum(double(cat(3,fr.cdata)),3)/n;
    I(I > 255) = 255;

    fr2.cdata = uint8(I);
    mov = addframe(mov,fr2);
    drawnow;
    
    timedWaitBar(ii/length(frames));
end;

mov = close(mov);
timedWaitBar(1);

