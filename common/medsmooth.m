function y = medsmooth(x, n)

sz = size(x);

if (mod(n,2) == 0),
  n = n+1;
end;

wdim = ndims(x)+1;
wnd = shiftdim((1:n)',-(wdim-1));
wnd = round(wnd - mean(wnd));
wnd = repmat(wnd,[sz 1]);

for dim = 1:wdim-1,
  i = (1:size(x,dim))';
  i = shiftdim(i,-(dim-1));

  sz1 = sz;
  sz1(dim) = 1;
  
  sub{dim} = repmat(i,[sz1 n]);
end;

sub{1} = sub{1} + wnd;
p = find((sub{1} < 1) | (sub{1} > size(x,1)));
sub{1}(p) = 1;

ind = sub2ind(size(x),sub{:});

val = x(ind);
val(p) = NaN;

y = nanmedian(val,wdim);

