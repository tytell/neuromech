%generate test matrices
sz = [100 100];
off = [3 2];
noisefrac = 0.1;

I1 = rand(sz);
I2 = rand(sz);

I2(off(1)+1:end,off(2)+1:end) = I1(1:end-off(1),1:end-off(2));

ind = round(rand(ceil(noisefrac*numel(I1)),1) * (numel(I1)-1)) + 1;

I2(ind) = I2(ind) + rand(size(ind));
I2(I2 > 1) = 1;

[x,y,u,v] = cpiv(I1,I2, [10 10 20 20 10 10]);



