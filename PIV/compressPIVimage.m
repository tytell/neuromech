function compressPIVimage(in,out, ratio)

nfrcalc = 100;

info = aviinfo(in);

N = info.NumFrames;

cmap = repmat((0:255)'/255,[1 3]);
i = 2:info.Height-1;
j = 2:info.Width-1;
npix = info.Height*info.Width;

timedWaitBar(0,'Compressing...');

mov = avifile(out,'compression','none','quality',100);
for fr = 1:N,
  I = frame2im(aviread(in,fr));

  Igrad = zeros(size(I));

  Id = im2double(I);
  Igrad(i,j) = (Id(i+1,j)-Id(i-1,j))/2 + ...
      (Id(i,j+1)-Id(i,j-1))/2;

  if (mod(fr-1,nfrcalc) == 0),
    p = prctile(Igrad(:), 100*(1-ratio));
  end;

  Ic = repmat(uint8(0),size(I));
  k = find(Igrad >= p);
  Ic(k) = I(k);

  actratio(fr) = length(k)/npix;

  F = im2frame(Ic,cmap);
  mov = addframe(mov,F);

  timedWaitBar(fr/N);
end;

mov = close(mov);

timedWaitBar(1);

fprintf('Mean compression ratio: %f\n', mean(actratio));
